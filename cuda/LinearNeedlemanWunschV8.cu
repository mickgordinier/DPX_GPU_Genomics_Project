#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"
#include "../c++/timing.h"

// Blocks are 1D with a size of the 32 threads (For 1 warp)
#define BLOCK_SIZE 32
#define BATCH_SIZE 100

// Defining this will test all of the sequences in the input file
#define TEST_ALL

/*
    THINGS TO CONSIDER FOR OPTIMIZATION
    1. Complete removal of the scoring matrix altogether (Use of warp shuffling and shared memory)
    2. Using 16x2 DPX instructions to have a thread work on 2 cells concurrently

*/

// NEEDLEMAN WUNSCH BASELINE KERNEL
// SUPPORTS KERNELS WITH THREADS LESS THAN QUERY LENGTH
// KERNEL IS LAUNCHED WITH MULTIPLE BLOCKS
// EACH BLOCK HANDLES A SEQUENCE PAIR

__global__ void 
needleman_wunsch_kernel(
    int *similarityScores,
    directionMain **batchBacktrackMatrices,
    char *backtrackStringsRet, 
    const char *allSequences, const seqPair *allSequenceInfo,
    const int matchWeight, const int mismatchWeight, const int gapWeight,
    const int startingSequenceIdx, const int stringLengthMax)
{

    const int tid = threadIdx.x;

    extern __shared__ int warpEdgeScore[]; 

    // We are launching multiple blocks, each of a warp of threads
    // Each block handles their own sequence alignment
    // We index into the array to obtain the strings and length
    
    directionMain *backtrackMatrix = batchBacktrackMatrices[blockIdx.x];

    const int sequenceIdx = startingSequenceIdx + blockIdx.x;
    const seqPair sequenceInfo = allSequenceInfo[sequenceIdx];
    
    const char *queryString = allSequences + sequenceInfo.queryIdx;
    const char *referenceString = allSequences + sequenceInfo.referenceIdx;

    // The matrices are of size (queryLength + 1) * (referenceLength + 1)
    const int numRows = sequenceInfo.querySize + 1;
    const int numCols = sequenceInfo.referenceSize + 1;

    /* --- (BEGIN) INITIALIZING THE SCORING MATRIX --- */

    // Used for when a thread has to iterate over more than one col/row
    int elementIdx;

    // Initialize the top row
    // Writing in DRAM burst for faster updating
    elementIdx = tid;
    while(elementIdx < numCols) {
        backtrackMatrix[elementIdx] = QUERY_INSERTION;
        elementIdx += BLOCK_SIZE;
    }

    // Initialize the left col
    // NOT Writing in DRAM burst (slower)
    elementIdx = tid;
    while(elementIdx < numRows) {
        backtrackMatrix[elementIdx*numCols] = QUERY_DELETION;
        elementIdx += BLOCK_SIZE;
    }

    if (tid == 0) {
        backtrackMatrix[0] = NONE_MAIN;
    }

    // Need to ensure that all threads in the block complete filling up all the edges
    // Do not need to do syncthreads across each loop iteration as there is no dependencies
    __syncthreads();

    /* --- (END) INITIALIZING THE SCORING MATRIX --- */

    /* --- (BEGIN) POPULATING THE SCORING MATRIX -- */

    /*
    relative cell indices
    [00][01]
    [10][11]
    */

    int leftDiag = gapWeight*tid;
    int left = gapWeight*(tid+1);
    int up = gapWeight*(tid+1); 
    int tmp_left = left;

    // Going through all of the rows each thread has to do
    for (int stripeStart = 1; stripeStart < numRows; stripeStart+=BLOCK_SIZE){

        int row = stripeStart + tid;
        int largestScore;

        /* threads outside of bounds should abort */
        if (row >= numRows) return;

        leftDiag = gapWeight*(row - 1);
        left = gapWeight*(row);

        for (int col = 1; col < (numCols+numRows); ++col){
            
            int adj_col = col - tid;

            if (row == 1){
                leftDiag = gapWeight*(adj_col - 1);
                up = gapWeight*(adj_col);
            }

            /* for all but the first stripe, t0 must grab its diagonal and upper values from t31 */
            if (stripeStart > 1 && tid == 0 && adj_col < numCols){
                up = warpEdgeScore[adj_col];
                leftDiag = (adj_col == 1) ? gapWeight*(row - 1) : warpEdgeScore[adj_col - 1];
            }

            if (adj_col > 0 && adj_col < numCols){
                largestScore = 0;
                char queryChar = queryString[row-1];
                char referenceChar = referenceString[adj_col-1];
    
                directionMain cornerDirection = NONE_MAIN;
                bool pred;
                bool isMatch = (queryChar == referenceChar);
                cornerDirection = isMatch ? MATCH : MISMATCH;
    
                int matchMismatchScore = isMatch ? leftDiag + matchWeight : leftDiag + mismatchWeight;
                int queryDeletionScore = up + gapWeight;
                int queryInsertionScore = left + gapWeight;
    
                largestScore = __vibmax_s32(queryDeletionScore, matchMismatchScore, &pred);
                if (pred) cornerDirection = QUERY_DELETION;
                        
                largestScore = __vibmax_s32(queryInsertionScore, largestScore, &pred);
                if (pred) cornerDirection = QUERY_INSERTION;

                // scoringMatrix[row * numCols + adj_col] = largestScore;
                backtrackMatrix[row * numCols + adj_col] = cornerDirection;

                tmp_left = left;
                left = largestScore;

                /* last thread in warp stores its scores in shared memory for t0 to access */
                if (tid == 31){
                    warpEdgeScore[adj_col] = largestScore;
                }
            }

            /*  top value for thread n + 1 is thread n's largestScore (just calculated value)*/
            up = __shfl_up_sync(0xffffffff, largestScore, 1);

            /* left diag value for thread n + 1 is thread n's left value (previously calculated value) */
            leftDiag = __shfl_up_sync(0xffffffff, tmp_left, 1);
        }

        if (row == numRows-1) {
            similarityScores[blockIdx.x] = largestScore;
        }
    }

    /* --- (END) POPULATING THE SCORING MATRIX -- */

    /* --- (BEGIN) DETERMINING BACKTRACKING -- */

    // Starting at the end
    if (tid == 0) {

        int referenceStrIdx = (stringLengthMax * 3) * blockIdx.x + (stringLengthMax-1);
        int alignmentStrIdx = referenceStrIdx + stringLengthMax;
        int queryStrIdx = alignmentStrIdx + stringLengthMax;

        backtrackStringsRet[referenceStrIdx] = '\0';
        backtrackStringsRet[alignmentStrIdx] = '\0';
        backtrackStringsRet[queryStrIdx] = '\0';

        int currentMemoRow = numRows - 1;
        int currentMemoCol = numCols - 1;

        while ((currentMemoRow != 0) || (currentMemoCol != 0)) {

            referenceStrIdx--;
            alignmentStrIdx--;
            queryStrIdx--;
            
            // Determine the current cell's predecessor
            switch (backtrackMatrix[(currentMemoRow * numCols) + currentMemoCol]) {
                
                case MATCH:
                    backtrackStringsRet[referenceStrIdx] = referenceString[currentMemoCol-1];
                    backtrackStringsRet[alignmentStrIdx] = '*';
                    backtrackStringsRet[queryStrIdx] = queryString[currentMemoRow-1];
                    --currentMemoRow;
                    --currentMemoCol;
                    break;
                // end if match

                case MISMATCH: 
                    backtrackStringsRet[referenceStrIdx] = referenceString[currentMemoCol-1];
                    backtrackStringsRet[alignmentStrIdx] = '|';
                    backtrackStringsRet[queryStrIdx] = queryString[currentMemoRow-1];
                    --currentMemoRow;
                    --currentMemoCol;
                    break;
                // end if mismatch
                
                case QUERY_DELETION:
                    backtrackStringsRet[referenceStrIdx] = '_';
                    backtrackStringsRet[alignmentStrIdx] = ' ';
                    backtrackStringsRet[queryStrIdx] = queryString[currentMemoRow-1];
                    --currentMemoRow;
                    break;
                // end if query deletion
                
                case QUERY_INSERTION:
                    backtrackStringsRet[referenceStrIdx] = referenceString[currentMemoCol-1];
                    backtrackStringsRet[alignmentStrIdx] = ' ';
                    backtrackStringsRet[queryStrIdx] = '_';
                    --currentMemoCol;
                    break;
                // end if query insertion
                
                default:
                    printf("ERROR\n");
                    return;
                // end if upper gap

            } // end switch
        } // end while

        // printf("TEMP Sequence Idx: %d, Reference: %s\n", sequenceIdx, backtrackStringsRet + referenceStrIdx);
        // printf("TEMP Sequence Idx: %d, Alignment: %s\n", sequenceIdx, backtrackStringsRet + alignmentStrIdx);
        // printf("TEMP Sequence Idx: %d, Query    : %s\n", sequenceIdx, backtrackStringsRet + queryStrIdx);

        int spacing = referenceStrIdx - ((stringLengthMax * 3) * blockIdx.x);

        for (int i = 0; i < spacing; ++i) {
            backtrackStringsRet[--referenceStrIdx] = '\0';
            backtrackStringsRet[--alignmentStrIdx] = '\0';
            backtrackStringsRet[--queryStrIdx] = '\0';
        }
    }
}


void
handleErrs(
    cudaError_t err,
    const char *errMsg) 
{
    if (err != cudaSuccess) {
        printf(errMsg);
        printf("CUDA ERROR: %s\n", cudaGetErrorString(err));
        exit(1);
    }
}


int main(int argc, char *argv[]) {

    // Print some cuda details
    printf("[Cuda Details]\n");
    int deviceCount;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess) {
        printf("FAILED TO GET DEVICE COUNT\n");
        printf("CUDA test kernel error: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    printf("Device count: %d\n", deviceCount);
    int device = 0;
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    printf("Device %d has compute capability %d.%d.\n",
           device, deviceProp.major, deviceProp.minor);
    printf("Concurrent kernels?: %d\n\n", deviceProp.concurrentKernels);

    // Check that YOU use it correctly
    if (argc < 2) {
		fprintf(stderr, "usage: main -pairs <InSeqFile> -match <matchWeight> -mismatch <mismatchWeight> -gap <gapWeight> \n");
		exit(EXIT_FAILURE);
    }
	
    // Get args
    char *pairFileName;
    int matchWeight     = 3;
    int mismatchWeight  = -1;
    int gapWeight       = -2;
    if(strcmp(argv[1], "-pairs") == 0) {
        pairFileName = argv[2];
    }
    if(argc > 3 && strcmp(argv[3], "-match") == 0) {
        matchWeight = atoi(argv[4]);
    }
    if(argc > 5 && strcmp(argv[5], "-mismatch") == 0) {
        mismatchWeight = atoi(argv[6]);
    }
    if(argc > 7 && strcmp(argv[7], "-gap") == 0) {
        gapWeight = atoi(argv[8]);
    }

    // Parse input file
    printf("Parsing input file: %s\n", pairFileName);
    inputInfo fileInfo;
    seqPair* allSequenceInfo;
    char* sequences;
    fileInfo = parseInput(pairFileName, allSequenceInfo, sequences);
    printf("Num Pairs: %d\n\n", fileInfo.numPairs);

    // Start timer
    uint64_t kernel_time = 0;
    uint64_t memalloc_time = 0;
    uint64_t backtracking_time = 0;
    uint64_t start_time = start_timer();
    #ifdef TEST_ALL
        
        // Copy over the sequences
        char* deviceSequences;
        seqPair *deviceAllSequenceInfo;

        uint64_t start_memalloc = get_time();
        handleErrs(
            cudaMalloc(&deviceSequences, (fileInfo.numBytes) * sizeof(char)),
            "FAILED TO ALLOCATE MEMORY FOR ALL SEQUENCES\n"
        );

        handleErrs(
            cudaMemcpy(deviceSequences, sequences, (fileInfo.numBytes) * sizeof(char), cudaMemcpyHostToDevice),
            "FAILED TO COPY MEMORY FOR ALL SEQUENCES\n"
        );

        handleErrs(
            cudaMalloc(&deviceAllSequenceInfo, (fileInfo.numPairs) * sizeof(seqPair)),
            "FAILED TO ALLOCATE MEMORY FOR ALL SEQUENCES\n"
        );

        handleErrs(
            cudaMemcpy(deviceAllSequenceInfo, allSequenceInfo, (fileInfo.numPairs) * sizeof(seqPair), cudaMemcpyHostToDevice),
            "FAILED TO COPY MEMORY FOR ALL SEQUENCES\n"
        );

        directionMain **tempDeviceBacktrackMatrices = (directionMain **)malloc(BATCH_SIZE * sizeof(directionMain*));

        int *deviceSimilarityScores;
        directionMain **deviceBacktrackMatrices;

        int *hostSimilarityScores = new int[BATCH_SIZE];

        handleErrs(
            cudaMalloc(&deviceSimilarityScores, BATCH_SIZE * sizeof(int)),
            "FAILED TO ALLOCATE MEMORY TO deviceBacktrackMatrices\n"
        );

        handleErrs(
            cudaMalloc(&deviceBacktrackMatrices, BATCH_SIZE * sizeof(directionMain*)),
            "FAILED TO ALLOCATE MEMORY TO deviceBacktrackMatrices\n"
        );
        memalloc_time += get_time() - start_memalloc;

        // Run the kernel on every sequence
        for(size_t sequenceIdx = 0; sequenceIdx < 2000; sequenceIdx+=BATCH_SIZE){
            
            start_memalloc = get_time();

            int largestReferenceLength = 0;
            int largestQueryLength = 0;

            for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; ++i) {
                
                const int queryLength = allSequenceInfo[i].querySize;
                const int referenceLength = allSequenceInfo[i].referenceSize;

                largestReferenceLength = max(largestReferenceLength, referenceLength);
                largestQueryLength = max(largestQueryLength, referenceLength);
                
                directionMain *deviceBacktrackMatrix;
        
                handleErrs(
                    cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(directionMain)),
                    "FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n"
                );

                tempDeviceBacktrackMatrices[i-sequenceIdx] = deviceBacktrackMatrix;
            }

            handleErrs(
                cudaMemcpy(deviceBacktrackMatrices, tempDeviceBacktrackMatrices, BATCH_SIZE * sizeof(directionMain*), cudaMemcpyHostToDevice),
                "FAILED TO COPY MEMORY FOR deviceBacktrackMatrices\n"
            );

            int stringLengthMax = (largestReferenceLength+largestQueryLength+1);

            char *deviceBacktrackStringRet;
        
            handleErrs(
                cudaMalloc(&deviceBacktrackStringRet, (stringLengthMax * 3) * BATCH_SIZE * sizeof(char)),
                "FAILED TO ALLOCATE MEMORY TO BACKTRACKING STRINGS\n"
            );

            memalloc_time += get_time() - start_memalloc;

            uint64_t start_kernel = get_time();
            // Need to launch kernel
            int smem_size = (largestReferenceLength + 1) * sizeof(int);
            needleman_wunsch_kernel<<<BATCH_SIZE, BLOCK_SIZE, smem_size>>>(
                deviceSimilarityScores,
                deviceBacktrackMatrices,
                deviceBacktrackStringRet,
                deviceSequences, deviceAllSequenceInfo,
                matchWeight, mismatchWeight, gapWeight,
                sequenceIdx, stringLengthMax
            );
            
            // Wait for kernel to finish
            handleErrs(
                cudaDeviceSynchronize(),
                "SYNCHRONIZATION FAILED\n"
            );
            kernel_time += get_time() - start_kernel;

            start_memalloc = get_time();

            handleErrs(
                cudaMemcpy(hostSimilarityScores, deviceSimilarityScores, BATCH_SIZE * sizeof(int), cudaMemcpyDeviceToHost),
                "FAILED TO COPY SIMILARITY SCORES FROM DEVICE --> HOST\n"
            );

            char *hostBacktrackingStringRet = (char *)malloc(stringLengthMax * (BATCH_SIZE) * sizeof(char));

            handleErrs(
                cudaMemcpy(hostBacktrackingStringRet, deviceBacktrackStringRet, (stringLengthMax * 3) * BATCH_SIZE * sizeof(char), cudaMemcpyDeviceToHost),
                "FAILED TO COPY SIMILARITY SCORES FROM DEVICE --> HOST\n"
            );

            memalloc_time += get_time() - start_memalloc;

            for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; ++i) {
                
                // start_memalloc = get_time();
                // const char *queryString = &sequences[allSequenceInfo[i].queryIdx];
                // const char *referenceString = &sequences[allSequenceInfo[i].referenceIdx];

                // const int queryLength = allSequenceInfo[i].querySize;
                // const int referenceLength = allSequenceInfo[i].referenceSize;

                // Copy the matrices back over
                // directionMain *hostBacktrackMatrix = new directionMain[(referenceLength+1) * (queryLength+1)];

                // Copy information back from device --> host
                // handleErrs(
                //     cudaMemcpy(hostBacktrackMatrix, tempDeviceBacktrackMatrices[i-sequenceIdx], (referenceLength+1) * (queryLength+1) * sizeof(directionMain), cudaMemcpyDeviceToHost),
                //     "FAILED TO COPY BACKTRACK MATRIX FROM DEVICE --> HOST\n"
                // );

                // // cudaFree(tempDeviceScoringMatrices[i-sequenceIdx]);
                // cudaFree(tempDeviceBacktrackMatrices[i-sequenceIdx]);
                // memalloc_time += get_time() - start_memalloc;


                // Backtrack matrices
                printf("%d | %d\n", i, hostSimilarityScores[i-sequenceIdx]);
                // uint64_t start_backtrack = get_time();
                // backtrackNW(hostBacktrackMatrix, referenceString, referenceLength, queryString, queryLength);
                // backtracking_time += get_time() - start_backtrack;

                // Free data arrays
                // delete[] hostScoringMatrix;
                // delete[] hostBacktrackMatrix;

                int spacing = ((stringLengthMax * 3) * (i-sequenceIdx));

                while (hostBacktrackingStringRet[((stringLengthMax * 3) * (i-sequenceIdx)) + spacing] == '\0') {
                    ++spacing;
                }

                printf("Sequence Idx: %d, Reference: %s\n", i, hostBacktrackingStringRet + ((stringLengthMax * 3) * (i-sequenceIdx)) + spacing);
                printf("Sequence Idx: %d, Alignment: %s\n", i, hostBacktrackingStringRet + ((stringLengthMax * 3) * (i-sequenceIdx)) + stringLengthMax + spacing);
                printf("Sequence Idx: %d, Query    : %s\n", i, hostBacktrackingStringRet + ((stringLengthMax * 3) * (i-sequenceIdx)) + stringLengthMax + stringLengthMax + spacing);
            }

            free(hostBacktrackingStringRet);
        }

        cudaFree(deviceSequences);
        cudaFree(deviceAllSequenceInfo);

        // free(tempDeviceScoringMatrices);
        free(tempDeviceBacktrackMatrices);
        free(hostSimilarityScores);

        // cudaFree(deviceScoringMatrices);
        cudaFree(deviceBacktrackMatrices);
        cudaFree(deviceSimilarityScores);
    #else

        // Copy over the sequences
        char* deviceSequences;

        handleErrs(
            cudaMalloc(&deviceSequences, (fileInfo.numBytes) * sizeof(char)),
            "FAILED TO ALLOCATE MEMORY FOR ALL SEQUENCES\n"
        );

        handleErrs(
            cudaMemcpy(deviceSequences, sequences, (fileInfo.numBytes) * sizeof(char), cudaMemcpyHostToDevice),
            "FAILED TO COPY MEMORY FOR ALL SEQUENCES\n"
        );


        

    
        int *scoringMatrix,          directionMain *backtrackMatrix,
        const char *queryStrings,    const char *referenceStrings,
        const int *queryStartingIdxs, const int *referenceStartingIdxs,
        const int *queryLengths,     const int *referenceLengths,
        const int matchWeight,       const int mismatchWeight, 
        const int gapWeight

        char *











        
        
    char *referenceString = &sequences[sequenceIdxs[0].referenceIdx];
        char *queryString = &sequences[sequenceIdxs[0].queryIdx];
        // char *referenceString = "GTCATGCAATAACG";
        // char *queryString = "ATGCAATA";
        // char *referenceString = "GTCAGTA";
        // char *queryString = "ATACA";

        int referenceLength = strlen(referenceString);
        int queryLength = strlen(queryString);

        printf("Reference String: %s (Length: %d)\n", referenceString, referenceLength);
        printf("Query String: %s (Length: %d)\n", queryString, queryLength);
        printf("(MATCH WEIGHT, MISMATCH WEIGHT, GAP WEIGHT): (%d, %d, %d)\n\n", matchWeight, mismatchWeight, gapWeight);

        // Allocate device memory for matrices
        printf("[Allocating CUDA Memory]\n");
        int *deviceScoringMatrix;
        directionMain *deviceBacktrackMatrix;
        char *deviceReferenceString;
        char *deviceQueryString;

        handleErrs(
            cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)),
            "FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n"
        );

        handleErrs(
            cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(directionMain)),
            "FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n"
        );

        handleErrs(
            cudaMalloc(&deviceReferenceString, (referenceLength) * sizeof(char)),
            "FAILED TO ALLOCATE MEMORY TO REFERENCE STRING\n"
        );

        handleErrs(
            cudaMemcpy(deviceReferenceString, referenceString, (referenceLength) * sizeof(char), cudaMemcpyHostToDevice),
            "FAILED TO COPY MEMORY TO REFERENCE STRING\n"
        );

        handleErrs(
            cudaMalloc(&deviceQueryString, (queryLength) * sizeof(char)),
            "FAILED TO ALLOCATE MEMORY TO QUERY STRING\n"
        );

        handleErrs(
            cudaMemcpy(deviceQueryString, queryString, (queryLength) * sizeof(char), cudaMemcpyHostToDevice),
            "FAILED TO COPY MEMORY TO QUERY STRING\n"
        );

        // Need to launch sinular kernel
        // Launching a kernel with 1 block with BLOCK_SIZE threads to populate scoring matrix
        needleman_wunsch_kernel<<<1, BLOCK_SIZE>>>(
            deviceScoringMatrix, deviceBacktrackMatrix,
            deviceQueryString, deviceReferenceString, 
            queryLength, referenceLength, 
            matchWeight, mismatchWeight, gapWeight
        );

        // Wait for kernel to finish
        handleErrs(
            cudaDeviceSynchronize(),
            "SYNCHRONIZATION FAILED\n"
        );

        // Allocate host memory for matrices
        // Allow for matrices to come from device -> host
        // Free up device memory
        int *hostScoringMatrix = new int[(referenceLength+1) * (queryLength+1)];
        directionMain *hostBacktrackMatrix = new directionMain[(referenceLength+1) * (queryLength+1)];

        // Copy information back from device --> host
        handleErrs(
            cudaMemcpy(hostScoringMatrix, deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost),
            "FAILED TO COPY SCORING MATRIX FROM DEVICE --> HOST"
        );
        
        handleErrs(
            cudaMemcpy(hostBacktrackMatrix, deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(directionMain), cudaMemcpyDeviceToHost),
            "FAILED TO COPY BACKTRACK MATRIX FROM DEVICE --> HOST"
        );

        cudaFree(deviceScoringMatrix);
        cudaFree(deviceBacktrackMatrix);
        cudaFree(deviceQueryString);
        cudaFree(deviceReferenceString);

        // Print Matrix
        printf("Scored Matrix\n");
        printMatrix(hostScoringMatrix, referenceLength + 1, queryLength + 1);
        printf("Backtrack Matrix\n");
        printBacktrackMatrix(hostBacktrackMatrix, referenceLength + 1, queryLength + 1);
        

        // Perform backtracking on host and print results
        printf("0 | %d\n", hostScoringMatrix[(referenceLength + 1) * (queryLength + 1) - 1]);
        uint64_t start_backtrack = get_time();
        backtrackNW(hostBacktrackMatrix, referenceString, referenceLength, queryString, queryLength);
        backtracking_time += get_time() - start_backtrack;
        
        // Free data arrays
        delete[] hostScoringMatrix;
        delete[] hostBacktrackMatrix;
    #endif

    uint64_t elapsed_time = get_elapsed_time();
    printf("Elapsed time (usec): %lld\n", elapsed_time);
    printf("Elapsed kernel time (usec): %lld\n", kernel_time);
    printf("Elapsed backtracking time (usec): %lld\n", backtracking_time);
    printf("Elapsed memallocing time (usec): %lld\n", memalloc_time);
    printf("Elapsed time sum (usec): %lld\n",kernel_time + backtracking_time + memalloc_time);

    // Cleanup
    printf("Cleaning up\n");
    cleanupParsedFile(allSequenceInfo, sequences);
}