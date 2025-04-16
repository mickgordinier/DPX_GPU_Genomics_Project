#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"

// Blocks are 1D with a size of the 32 threads (For 1 warp)
#define BLOCK_SIZE 32

// Defing this will test all of the sequences in the input file
#define TEST_ALL

__global__ void 
needleman_wunsch_kernel_warp_shuffle(
    int * scoringMatrix,
    direction * __restrict__ backtrackMatrix,
    const char * __restrict__ queryString, const char * __restrict__ referenceString,
    const int queryLength, const int referenceLength,
    const int matchWeight, const int mismatchWeight, 
    const int gapWeight)
{
    // We are only launching 1 block
    // Thus, each thread will only have a unique threadID that differentiates the threads

    const int tid = threadIdx.x;
    const int threadCount = blockDim.x;

    const int numRows = queryLength + 1;
    const int numCols = referenceLength + 1;

    extern __shared__ int warpEdgeScore[];

    /* --- (BEGIN) INITIALIZING THE SCORING MATRIX --- */
    // Used for when a thread has to iterate over more than one col/row
    int elementIdx = tid;

    // Initialize the top row
    // Writing in DRAM burst for faster updating

    while(elementIdx < numCols) {
        scoringMatrix[elementIdx] = gapWeight*elementIdx;
        backtrackMatrix[elementIdx] = QUERY_INSERTION;
        elementIdx += threadCount;
    }

    // Initialize the left col
    // NOT Writing in DRAM burst (slower)
    elementIdx = tid;
    while(elementIdx < numRows) {
        scoringMatrix[elementIdx*numCols] = gapWeight*elementIdx;
        backtrackMatrix[elementIdx*numCols] = QUERY_DELETION;
        elementIdx += threadCount;
    }

    if (tid == 0) {
        backtrackMatrix[0] = NONE;
    }

    /* --- (END) INITIALIZING THE SCORING MATRIX --- */

    /* --- (BEGIN) POPULATING THE SCORING MATRIX -- */
    int leftDiag = gapWeight*tid, left = gapWeight*(tid+1), up = gapWeight*(tid+1), tmp_left = left;

    for (int stripeStart = 1; stripeStart < numRows; stripeStart+=threadCount){

        int row = stripeStart + tid;

        /* threads outside of bounds should abort */
        if (row >= numRows) return;

        for (int col = 1; col < (numCols+numRows); ++col){
            int largestScore = 0;
            int adj_col = col - tid;

            if (row == 1){
                leftDiag = gapWeight*(adj_col - 1);
                up = gapWeight*(adj_col);
            }

            if (adj_col == 1){
                leftDiag = gapWeight*(row - 1);
                left = gapWeight*(row);
            }

            /* for all but the first stripe, t0 must grab its diagonal and upper values from t31 */
            if (stripeStart > 1 && tid == 0 && adj_col < numCols){
                up = warpEdgeScore[adj_col];
                leftDiag = (adj_col == 1) ? gapWeight*(row - 1) : warpEdgeScore[adj_col - 1];
            }

            if (adj_col > 0 && adj_col < numCols){
                char queryChar = queryString[row-1];
                char referenceChar = referenceString[adj_col-1];
    
                direction cornerDirection = NONE;
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

                scoringMatrix[row * numCols + adj_col] = largestScore;
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
    int threadCount     = 32;
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
    if(argc > 9 && strcmp(argv[9], "-threads-per-alignment") == 0) {
        threadCount = atoi(argv[10]);
    }

    // Parse input file
    printf("Parsing input file: %s\n", pairFileName);
    inputInfo fileInfo;
    seqPair* sequenceIdxs;
    char* sequences;
    fileInfo = parseInput(pairFileName, sequenceIdxs, sequences);
    printf("Num Pairs: %d\n\n", fileInfo.numPairs);

    #ifdef TEST_ALL
        
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

        // Run the kernel on every sequence
        for(size_t i = 0; i < fileInfo.numPairs; i++){

            char *referenceString = &sequences[sequenceIdxs[i].referenceIdx];
            char *queryString = &sequences[sequenceIdxs[i].queryIdx];

            int referenceLength = strlen(referenceString);
            int queryLength = strlen(queryString);

            int *deviceScoringMatrix;
            direction *deviceBacktrackMatrix;

            handleErrs(
                cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)),
                "FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n"
            );
    
            handleErrs(
                cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction)),
                "FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n"
            );

            // Need to launch kernel
            /* allocate enough shared memory to store 1 row of scores */
            int smem_size = (referenceLength+1)*sizeof(int);
            needleman_wunsch_kernel_warp_shuffle<<<1, BLOCK_SIZE, smem_size>>>(
                deviceScoringMatrix, deviceBacktrackMatrix,
                deviceSequences + sequenceIdxs[i].queryIdx, deviceSequences + sequenceIdxs[i].referenceIdx, 
                sequenceIdxs[i].querySize, sequenceIdxs[i].referenceSize, 
                matchWeight, mismatchWeight, gapWeight
            );
            
            // Wait for kernel to finish
            handleErrs(
                cudaDeviceSynchronize(),
                "SYNCHRONIZATION FAILED\n"
            );

            // Copy the matrices back over
            int *hostScoringMatrix = new int[(referenceLength+1) * (queryLength+1)];
            direction *hostBacktrackMatrix = new direction[(referenceLength+1) * (queryLength+1)];

            // Copy information back from device --> host
            handleErrs(
                cudaMemcpy(hostScoringMatrix, deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost),
                "FAILED TO COPY SCORING MATRIX FROM DEVICE --> HOST"
            );
            
            handleErrs(
                cudaMemcpy(hostBacktrackMatrix, deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction), cudaMemcpyDeviceToHost),
                "FAILED TO COPY BACKTRACK MATRIX FROM DEVICE --> HOST"
            );

            cudaFree(deviceScoringMatrix);
            cudaFree(deviceBacktrackMatrix);

            // Backtrack matrices
            printf("%d | %d\n", i, hostScoringMatrix[(referenceLength + 1) * (queryLength + 1) - 1]);
            backtrackNW(hostBacktrackMatrix, referenceString, referenceLength, queryString, queryLength);

            // Free data arrays
            delete[] hostScoringMatrix;
            delete[] hostBacktrackMatrix;
        }

        cudaFree(deviceSequences);
    #else
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
        direction *deviceBacktrackMatrix;
        char *deviceReferenceString;
        char *deviceQueryString;

        handleErrs(
            cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)),
            "FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n"
        );

        handleErrs(
            cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction)),
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
        // Launching a kernel with 1 block with threadCount threads to populate scoring matrix
        int smem_size = (referenceLength+1)*sizeof(int);
        needleman_wunsch_kernel_warp_shuffle<<<1, BLOCK_SIZE, smem_size>>>(
            deviceScoringMatrix, deviceBacktrackMatrix,
            deviceSequences + sequenceIdxs[i].queryIdx, deviceSequences + sequenceIdxs[i].referenceIdx, 
            sequenceIdxs[i].querySize, sequenceIdxs[i].referenceSize, 
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
        direction *hostBacktrackMatrix = new direction[(referenceLength+1) * (queryLength+1)];

        // Copy information back from device --> host
        handleErrs(
            cudaMemcpy(hostScoringMatrix, deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost),
            "FAILED TO COPY SCORING MATRIX FROM DEVICE --> HOST"
        );
        
        handleErrs(
            cudaMemcpy(hostBacktrackMatrix, deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction), cudaMemcpyDeviceToHost),
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
        backtrackNW(hostBacktrackMatrix, referenceString, referenceLength, queryString, queryLength);
        
        // Free data arrays
        delete[] hostScoringMatrix;
        delete[] hostBacktrackMatrix;
    #endif

    // Cleanup
    printf("Cleaning up\n");
    cleanupParsedFile(sequenceIdxs, sequences);
}