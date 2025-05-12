#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"
#include "../c++/timing.h"

// Blocks are 1D with a size of the 32 threads (For 1 warp)
#define BLOCK_SIZE 32
#define BATCH_SIZE 1000


__global__ void 
needleman_wunsch_kernel(
    unsigned int *similarityScores,
    int *stringSpacing,
    directionMain *batchBacktrackMatrices,
    int *batchBacktrackMatrixStartingIndices,
    char *backtrackStringsRet, 
    const char *allSequences, const seqPair *allSequenceInfo,
    const unsigned int matchWeight, const unsigned int mismatchWeight, const unsigned int gapWeight,
    const int startingSequenceIdx, const int stringLengthMax)
{

    const int tid = threadIdx.x;

    extern __shared__ int warpEdgeScore[]; 

    // We are launching multiple blocks, each of a warp of threads
    // Each block handles their own sequence alignment
    // We index into the array to obtain the strings and length
    
    int startingMatrixIdx1 = batchBacktrackMatrixStartingIndices[(2 * blockIdx.x)];
    int startingMatrixIdx2 = batchBacktrackMatrixStartingIndices[(2 * blockIdx.x) + 1];

    directionMain *backtrackMatrix1 = batchBacktrackMatrices + startingMatrixIdx1;
    directionMain *backtrackMatrix2 = batchBacktrackMatrices + startingMatrixIdx2;

    const int sequenceIdx1 = startingSequenceIdx + (2 * blockIdx.x) ;
    const int sequenceIdx2 = startingSequenceIdx + (2 * blockIdx.x) + 1;

    const seqPair sequenceInfo1 = allSequenceInfo[sequenceIdx1];
    const seqPair sequenceInfo2 = allSequenceInfo[sequenceIdx2];

    const char *queryString1 = allSequences + sequenceInfo1.queryIdx;
    const char *queryString2 = allSequences + sequenceInfo2.queryIdx;

    const char *referenceString1 = allSequences + sequenceInfo1.referenceIdx;
    const char *referenceString2 = allSequences + sequenceInfo2.referenceIdx;

    const int numRows1 = sequenceInfo1.querySize + 1;
    const int numRows2 = sequenceInfo2.querySize + 1;

    const int numCols1 = sequenceInfo1.referenceSize + 1;
    const int numCols2 = sequenceInfo2.referenceSize + 1;

    /* --- (BEGIN) INITIALIZING THE SCORING MATRIX --- */

    // Used for when a thread has to iterate over more than one col/row
    int elementIdx;

    // Initialize the top row
    // Writing in DRAM burst for faster updating
    elementIdx = tid;
    while(elementIdx < numCols1) {
        backtrackMatrix1[elementIdx] = QUERY_INSERTION;
        elementIdx += BLOCK_SIZE;
    }
    
    elementIdx = tid;
    while(elementIdx < numCols2) {
        backtrackMatrix2[elementIdx] = QUERY_INSERTION;
        elementIdx += BLOCK_SIZE;
    }

    // Initialize the left col
    // NOT Writing in DRAM burst (slower)
    elementIdx = tid;
    while(elementIdx < numRows1) {
        backtrackMatrix1[elementIdx*numCols1] = QUERY_DELETION;
        elementIdx += BLOCK_SIZE;
    }

    elementIdx = tid;
    while(elementIdx < numRows2) {
        backtrackMatrix2[elementIdx*numCols2] = QUERY_DELETION;
        elementIdx += BLOCK_SIZE;
    }

    if (tid == 0) {
        backtrackMatrix1[0] = NONE_MAIN;
        backtrackMatrix2[0] = NONE_MAIN;
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

    unsigned int leftDiag = gapWeight * tid;
    unsigned int left     = gapWeight * (tid + 1);
    unsigned int up       = gapWeight * (tid + 1); 

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

                left = largestScore;

                /* last thread in warp stores its scores in shared memory for t0 to access */
                if (tid == 31){
                    warpEdgeScore[adj_col] = largestScore;
                }

                leftDiag = up;
            }

            /*  top value for thread n + 1 is thread n's largestScore (just calculated value)*/
            up = __shfl_up_sync(0xffffffff, largestScore, 1);
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

        stringSpacing[blockIdx.x] = referenceStrIdx;
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

    /* 
    store all backtracking matrices for a batch in one matrix - each warp will index in via index array
    */
    int *deviceBacktrackingIndices;
    int *hostBacktrackingIndices = (int *)malloc(BATCH_SIZE * sizeof(int));
    
    handleErrs(
        cudaMalloc(&deviceBacktrackingIndices, BATCH_SIZE * sizeof(int)),
        "FAILED TO ALLOCATE MEMORY TO deviceBacktrackingIndices\n"
    );
    
    int *deviceSimilarityScores;
    int *hostSimilarityScores = (int*)malloc(BATCH_SIZE * sizeof(int));

    handleErrs(
        cudaMalloc(&deviceSimilarityScores, BATCH_SIZE * sizeof(int)),
        "FAILED TO ALLOCATE MEMORY TO deviceSimilarityScores\n"
    );
    
    int *deviceStringSpacing;
    int *hostStringSpacing = (int*)malloc(BATCH_SIZE * sizeof(int));

    handleErrs(
        cudaMalloc(&deviceStringSpacing, BATCH_SIZE * sizeof(int)),
        "FAILED TO ALLOCATE MEMORY TO deviceSimilarityScores\n"
    );

    memalloc_time += get_time() - start_memalloc;

    // Run the kernel on every sequence
    for(size_t sequenceIdx = 0; sequenceIdx < fileInfo.numPairs; sequenceIdx+=BATCH_SIZE){
        start_memalloc = get_time();

        int largestReferenceLength = 0;
        int largestQueryLength = 0;

        /* first warp's starting index is 0 */
        uint64_t batchMatrixSize = 0;
        hostBacktrackingIndices[0] = 0;
        for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; ++i) {
            const int queryLength = allSequenceInfo[i].querySize;
            const int referenceLength = allSequenceInfo[i].referenceSize;

            largestReferenceLength = max(largestReferenceLength, referenceLength);
            largestQueryLength = max(largestQueryLength, queryLength);

            /* make sure we don't go over the end of the array */
            batchMatrixSize += ((referenceLength + 1) * (queryLength + 1));
            if ((i - sequenceIdx) < (BATCH_SIZE - 1)){
                hostBacktrackingIndices[i-sequenceIdx + 1] = batchMatrixSize;
            }
        }

        /* copy backtracking indices to device */
        handleErrs(
            cudaMemcpy(deviceBacktrackingIndices, hostBacktrackingIndices, BATCH_SIZE * sizeof(int), cudaMemcpyHostToDevice),
            "FAILED TO COPY MEMORY FOR deviceBacktrackingIndices\n"
        );

        /* allocate device mem for all backtracking matrices */
        directionMain *deviceMatricesAll;
        handleErrs(
            cudaMalloc(&deviceMatricesAll, batchMatrixSize*sizeof(directionMain)),
            "FAILED TO ALLOCATE MEMORY TO deviceMatricesAll\n"
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
            deviceStringSpacing,
            deviceMatricesAll,
            deviceBacktrackingIndices,
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

        handleErrs(
            cudaMemcpy(hostStringSpacing, deviceStringSpacing, BATCH_SIZE * sizeof(int), cudaMemcpyDeviceToHost),
            "FAILED TO COPY SIMILARITY SCORES FROM DEVICE --> HOST\n"
        );

        char *hostBacktrackingStringRet = (char *)malloc(stringLengthMax * 3 * BATCH_SIZE * sizeof(char));

        handleErrs(
            cudaMemcpy(hostBacktrackingStringRet, deviceBacktrackStringRet, (stringLengthMax * 3) * BATCH_SIZE * sizeof(char), cudaMemcpyDeviceToHost),
            "FAILED TO COPY BACKTRACKING STRING FROM DEVICE --> HOST\n"
        );

        memalloc_time += get_time() - start_memalloc;

        for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; ++i) {
        
            // Backtrack matrices
            printf("%d | %d\n", i, hostSimilarityScores[i-sequenceIdx]);

            int spacing = hostStringSpacing[i-sequenceIdx];

            printf("%s\n", hostBacktrackingStringRet + spacing);
            printf("%s\n", hostBacktrackingStringRet + stringLengthMax + spacing);
            printf("%s\n", hostBacktrackingStringRet + stringLengthMax + stringLengthMax + spacing);
        }

        free(hostBacktrackingStringRet);
        cudaFree(deviceBacktrackStringRet);
        cudaFree(deviceMatricesAll);
    }

    cudaFree(deviceSequences);
    cudaFree(deviceAllSequenceInfo);

    free(hostBacktrackingIndices);
    free(hostSimilarityScores);
    free(hostStringSpacing);

    cudaFree(deviceBacktrackingIndices);
    cudaFree(deviceSimilarityScores);
    cudaFree(deviceStringSpacing);


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