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
    int *similarityScores,
    int *stringSpacing,
    directionMain *batchBacktrackMatrices,
    int *batchIndices,
    int *stringStartingIndices,
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
    
    int batchIndex = batchIndices[blockIdx.x];
    directionMain *backtrackMatrix = batchBacktrackMatrices + batchIndex;

    const int sequenceIdx = startingSequenceIdx + blockIdx.x;
    const seqPair sequenceInfo = allSequenceInfo[sequenceIdx];
    
    const char *queryString = allSequences + sequenceInfo.queryIdx;
    const char *referenceString = allSequences + sequenceInfo.referenceIdx;

    const int queryLength = sequenceInfo.querySize;
    const int referenceLength = sequenceInfo.referenceSize;

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

        int jumpToNextString = queryLength + referenceLength + 1;

        int referenceStrIdx = stringStartingIndices[blockIdx.x] + (jumpToNextString - 1);
        int alignmentStrIdx = referenceStrIdx + jumpToNextString;
        int queryStrIdx = alignmentStrIdx + jumpToNextString;

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

    int *deviceStringStartingIndices;
    int *hostStringStartingIndices = (int *)malloc(BATCH_SIZE * sizeof(int));
    
    handleErrs(
        cudaMalloc(&deviceStringStartingIndices, BATCH_SIZE * sizeof(int)),
        "FAILED TO ALLOCATE MEMORY TO hostStringStartingIndices\n"
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
        uint64_t totalStringsSize = 0;
        hostBacktrackingIndices[0] = 0;
        hostStringStartingIndices[0] = 0;
        for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; ++i) {
            const int queryLength = allSequenceInfo[i].querySize;
            const int referenceLength = allSequenceInfo[i].referenceSize;

            largestReferenceLength = max(largestReferenceLength, referenceLength);
            largestQueryLength = max(largestQueryLength, queryLength);

            /* make sure we don't go over the end of the array */
            batchMatrixSize += ((referenceLength + 1) * (queryLength + 1));
            totalStringsSize += (3 * (queryLength + referenceLength + 1));

            if ((i - sequenceIdx) < (BATCH_SIZE - 1)){
                hostBacktrackingIndices[i-sequenceIdx + 1] = batchMatrixSize;
                hostStringStartingIndices[i-sequenceIdx + 1] = totalStringsSize;
            }
        }

        /* copy backtracking indices to device */
        handleErrs(
            cudaMemcpy(deviceBacktrackingIndices, hostBacktrackingIndices, BATCH_SIZE * sizeof(int), cudaMemcpyHostToDevice),
            "FAILED TO COPY MEMORY FOR deviceBacktrackingIndices\n"
        );

        handleErrs(
            cudaMemcpy(deviceStringStartingIndices, hostStringStartingIndices, BATCH_SIZE * sizeof(int), cudaMemcpyHostToDevice),
            "FAILED TO COPY MEMORY FOR deviceStringStartingIndices\n"
        );

        /* allocate device mem for all backtracking matrices */
        directionMain *deviceMatricesAll;
        handleErrs(
            cudaMalloc(&deviceMatricesAll, batchMatrixSize*sizeof(directionMain)),
            "FAILED TO ALLOCATE MEMORY TO deviceMatricesAll\n"
        );

        char *deviceBacktrackStringRet;
        handleErrs(
            cudaMalloc(&deviceBacktrackStringRet, totalStringsSize * sizeof(char)),
            "FAILED TO ALLOCATE MEMORY TO BACKTRACKING STRINGS\n"
        );

        int stringLengthMax = (largestReferenceLength+largestQueryLength+1);

        memalloc_time += get_time() - start_memalloc;

        uint64_t start_kernel = get_time();
        // Need to launch kernel
        int smem_size = (largestReferenceLength + 1) * sizeof(int);
        needleman_wunsch_kernel<<<BATCH_SIZE, BLOCK_SIZE, smem_size>>>(
            deviceSimilarityScores,
            deviceStringSpacing,
            deviceMatricesAll,
            deviceBacktrackingIndices,
            deviceStringStartingIndices,
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

        char *hostBacktrackingStringRet = (char *)malloc(totalStringsSize * sizeof(char));

        handleErrs(
            cudaMemcpy(hostBacktrackingStringRet, deviceBacktrackStringRet, totalStringsSize * sizeof(char), cudaMemcpyDeviceToHost),
            "FAILED TO COPY BACKTRACKING STRING FROM DEVICE --> HOST\n"
        );

        memalloc_time += get_time() - start_memalloc;

        for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; ++i) {
        
            // Backtrack matrices
            printf("%d | %d\n", i, hostSimilarityScores[i-sequenceIdx]);

            const int queryLength = allSequenceInfo[i].querySize;
            const int referenceLength = allSequenceInfo[i].referenceSize;

            const int jumpToNextString = queryLength + referenceLength + 1;

            const int spacing = hostStringSpacing[i-sequenceIdx];

            printf("%s\n", hostBacktrackingStringRet + spacing);
            printf("%s\n", hostBacktrackingStringRet + jumpToNextString + spacing);
            printf("%s\n", hostBacktrackingStringRet + jumpToNextString + jumpToNextString + spacing);
        }

        free(hostBacktrackingStringRet);
        cudaFree(deviceBacktrackStringRet);
        cudaFree(deviceMatricesAll);
    }

    cudaFree(deviceSequences);
    cudaFree(deviceAllSequenceInfo);

    free(hostBacktrackingIndices);
    free(hostStringStartingIndices);
    free(hostSimilarityScores);
    free(hostStringSpacing);

    cudaFree(deviceBacktrackingIndices);
    cudaFree(deviceStringStartingIndices);
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