#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"

// Blocks are 1D with a size of the maximum 1024 threads
#define BLOCK_SIZE 1024

// Defing this will test all of the sequences in the input file
#define TEST_ALL

/*
    THINGS TO CONSIDER FOR OPTIMIZATION
    1. Complete removal of the scoring matrix altogether (Use of warp shuffling and shared memory)
    2. Using 16x2 DPX instructions to have a thread work on 2 cells concurrently

*/

// NEEDLEMAN WUNSCH BASELINE KERNEL
__global__ void needleman_wunsch_kernel(int *scoringMatrix, direction *backtrackMatrix,
                        const char *queryString, const char *referenceString,
                        const int queryLength, const int referenceLength,
                        const int matchWeight, const int mismatchWeight, 
                        const int gapWeight,   const int threadCount){
    const int tid = threadIdx.x;
    const int numRows = queryLength + 1;
    const int numCols = referenceLength + 1;

    // Initialize the matrix
    if (tid < numRows*numCols) {
        scoringMatrix[tid] = -999; // sentry value
    }
    __syncthreads();

    // Adjusted tid, used for when a thread has to iterate over more than one col/row
    int adjtid;

    // Initialize the top row
    adjtid = tid;
    while(adjtid < numCols) {
        scoringMatrix[adjtid] = gapWeight*adjtid;
        adjtid += threadCount;
    }

    // Initialize the left col
    adjtid = tid;
    while(adjtid < numRows) {
        scoringMatrix[adjtid*numCols] = gapWeight*adjtid;
        adjtid += threadCount;
    }
    __syncthreads();

    /*
    relative cell indices
    [00][01]
    [10][11]
    */

    // Every thread gets a row and char
    int rowIdx = tid + 1;
    int colIdxAdj = 0; // Adjustment cal in case a thread needs to do multiple rows
    char queryChar = queryString[rowIdx - 1];

    // Main loop
    const int numLoops = (numCols + numRows) + (numCols) * (queryLength/threadCount);
    for(int idx = 1; idx < numLoops; idx++){
        // Make sure we dont use any threads that would go out of bounds
        if(rowIdx < numRows){
            int colIdx = idx - tid + colIdxAdj; // Stagger threads by their thread id
            int cell00Idx;
            int cell01Idx;
            int cell10Idx;
            int cell11Idx;

            // Setup cell indices once a thread can start executing
            if(colIdx == 1){
                cell00Idx = (rowIdx-1)*numCols + colIdx - 1;
                cell01Idx = (rowIdx-1)*numCols + colIdx;
                cell10Idx = rowIdx*numCols + colIdx - 1;
                cell11Idx = rowIdx*numCols + colIdx; 
            } 

            // Main cell updating
            if((colIdx > 0) && (colIdx < numCols)){
                char referenceChar = referenceString[colIdx - 1];
                direction cornerDirection = NONE;
                bool pred;
                
                // Determine if match
                bool isMatch = (queryChar == referenceChar);
                cornerDirection = isMatch ? MATCH : MISMATCH;

                // Get all the possible scores
                int matchMismatchScore = isMatch ? scoringMatrix[cell00Idx] + matchWeight : scoringMatrix[cell00Idx] + mismatchWeight;
                int queryDeletionScore = scoringMatrix[cell01Idx] + gapWeight;
                int queryInsertionScore = scoringMatrix[cell10Idx] + gapWeight;

                // Find the largest of the 3 scores
                int largestScore;
                largestScore = __vibmax_s32(queryDeletionScore, matchMismatchScore, &pred);
                if (pred) cornerDirection = QUERY_DELETION;
                
                largestScore = __vibmax_s32(queryInsertionScore, largestScore, &pred);
                if (pred) cornerDirection = QUERY_INSERTION;

                // Update scoring matrix
                scoringMatrix[cell11Idx] = largestScore;
                backtrackMatrix[cell11Idx] = cornerDirection;
                cell00Idx += 1;
                cell01Idx += 1;
                cell10Idx += 1;
                cell11Idx += 1;
            }

            // If thread already completed its row, ready it for another
            if(colIdx == (numCols - 1)){
                rowIdx += threadCount;
                colIdxAdj -= colIdx;
                queryChar = queryString[rowIdx - 1];
            }
        }
        __syncthreads();
    }
    return;

}

/*
// Device kernel that each thread will be executing to fill in its respective row
__global__ void
needleman_wunsch_forward_pass_kernel(int *scoringMatrix, direction *backtrackMatrix,
                        const char *queryString, const char *referenceString,
                        const int queryLength, const int referenceLength,
                        const int matchWeight, const int mismatchWeight, const int gapWeight)
{

    // Obtaining the 1D unique block and thread Id for the specific thread
    int tid = threadIdx.x;
    int numRows = queryLength + 1;
    int numCols = referenceLength + 1;

    // Initialize the top row
    if (tid < numCols) {
        scoringMatrix[tid] = gapWeight*tid;
    }
    __syncthreads();

    // Initialize the left col
    if (tid < numRows) {
        scoringMatrix[tid*numCols] = gapWeight*tid;
    }
    __syncthreads();

    int queryInsertionScore;
    int queryDeletionScore;
    int matchMismatchScore;
    int largestScore;

    direction cornerDirection;

    int whileCount = 0;

    int rowIdx, matrixIdx, rowAboveIdx;
    char queryChar;

    bool pred;

    while ((whileCount * BLOCK_SIZE) + 1) {

        rowIdx = (whileCount * BLOCK_SIZE) + tid + 1;

        matrixIdx = rowIdx * referenceLength;
        rowAboveIdx = ((rowIdx-1) * referenceLength) + 1;

        // Each thread initializing the 0th column of their row to 0
        if (rowIdx < queryLength) {
            scoringMatrix[matrixIdx++] = gapWeight * tid;
            queryChar = queryString[rowIdx - 1]; 
        }

        // Need to ensure all threads in the block have written to their respective location before continuing
        __syncthreads();

        // Need to fill up all the rows
        // Starting on column idx 1
        for (int i = 1; i < referenceLength + BLOCK_SIZE - 1; ++i) {

            // On initialization, the lower threads need to wait for the upper thread to begin
            // if i < tid, we can assume the thread has not yet calculated the above value
            
            // At the end, we don't want the thread to do any more computation at the end of its row
            // Thus, we will have the thread stop updating the matrix once i >= (referenceLength + tid)
            if ((rowIdx < queryLength) && (i >= tid) && (i < (referenceLength + tid))) {

                cornerDirection = NONE;

                const char referenceChar = referenceString[i - tid];

                // If match, add the match score from the corner
                if (queryChar == referenceChar){
                    matchMismatchScore = scoringMatrix[rowAboveIdx-1] + matchWeight;
                    cornerDirection = MATCH;
                }
                // Otherwise, add the penalty score from the corner
                else {
                    matchMismatchScore = scoringMatrix[rowAboveIdx-1] + mismatchWeight;
                    cornerDirection = MISMATCH;
                }
                
                // Calculate potential gap scores
                queryDeletionScore = scoringMatrix[rowAboveIdx] + gapWeight;
                queryInsertionScore = scoringMatrix[matrixIdx-1] + gapWeight;
                
                largestScore = __vibmax_s32(queryDeletionScore, matchMismatchScore, &pred);
                if (pred) cornerDirection = QUERY_DELETION;
                
                largestScore = __vibmax_s32(queryInsertionScore, largestScore, &pred);
                if (pred) cornerDirection = QUERY_INSERTION;

                scoringMatrix[matrixIdx] = largestScore;
                backtrackMatrix[matrixIdx] = cornerDirection;

                ++matrixIdx;
                ++rowAboveIdx;
            }

            // Need to ensure all threads in the block have written to their respective location before continuing
            __syncthreads();

        }

        ++whileCount;

    }
}
*/

int main(int argc, char *argv[]) {

    // Print some cuda details
    printf("[Cuda Details]\n");
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
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
    printf("Num Pairs: %d\n", fileInfo.numPairs);

    #ifdef TEST_ALL
        // Copy over the sequences
        char* deviceSequences;
        if (cudaMalloc(&deviceSequences, (fileInfo.numBytes) * sizeof(char)) != cudaSuccess) {
            printf("FAILED TO ALLOCATE MEMORY TO SCORING MATRIX FOR SEQUENCE SEQUENCES\n");
            exit(1);
        }

        cudaMemcpy(deviceSequences, sequences, (fileInfo.numBytes) * sizeof(char), cudaMemcpyHostToDevice);
        
        // Run the kernel on every sequence
        for(size_t i = 0; i < fileInfo.numPairs; i++){


            char *referenceString = &sequences[sequenceIdxs[i].referenceIdx];
            char *queryString = &sequences[sequenceIdxs[i].queryIdx];

            int referenceLength = strlen(referenceString);
            int queryLength = strlen(queryString);

            int *deviceScoringMatrix;
            direction *deviceBacktrackMatrix;

            if (cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)) != cudaSuccess) {
                printf("FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n");
                exit(1);
            }

            if (cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction)) != cudaSuccess) {
                printf("FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n");
                exit(1);
            }

            // Need to launch kernel
            needleman_wunsch_kernel<<<1, threadCount>>>(
                deviceScoringMatrix, deviceBacktrackMatrix,
                deviceSequences + sequenceIdxs[i].queryIdx, deviceSequences + sequenceIdxs[i].referenceIdx, 
                sequenceIdxs[i].querySize, sequenceIdxs[i].referenceSize, 
                matchWeight, mismatchWeight, gapWeight,
                threadCount);
            
            // Wait for kernel to finish
            cudaError_t err = cudaDeviceSynchronize();
            if (err != cudaSuccess) {
                printf("CUDA test kernel error: %s\n", cudaGetErrorString(err));
            }

            // Copy the matrices back over
            int *hostScoringMatrix = new int[(referenceLength+1) * (queryLength+1)];
            direction *hostBacktrackMatrix = new direction[(referenceLength+1) * (queryLength+1)];

            cudaMemcpy(hostScoringMatrix, deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(hostBacktrackMatrix, deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction), cudaMemcpyDeviceToHost);

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

        if (cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)) != cudaSuccess) {
            printf("FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n");
            exit(1);
        }

        if (cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction)) != cudaSuccess) {
            printf("FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n");
            exit(1);
        }

        if (cudaMalloc(&deviceReferenceString, (referenceLength) * sizeof(char)) != cudaSuccess) {
            printf("FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n");
            exit(1);
        }
        cudaMemcpy(deviceReferenceString, referenceString, (referenceLength) * sizeof(char), cudaMemcpyHostToDevice);

        if (cudaMalloc(&deviceQueryString, (queryLength) * sizeof(char)) != cudaSuccess) {
            printf("FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n");
            exit(1);
        }
        cudaMemcpy(deviceQueryString, queryString, (queryLength) * sizeof(char), cudaMemcpyHostToDevice);

        

        // Need to launch sinular kernel
        // Launching a kernel with 1 block with threadCount threads to populate scoring matrix
        needleman_wunsch_kernel<<<1, threadCount>>>(
            deviceScoringMatrix, deviceBacktrackMatrix,
            deviceQueryString, deviceReferenceString, 
            queryLength, referenceLength, 
            matchWeight, mismatchWeight, gapWeight,
            threadCount);


        // Wait for kernels to finish
        cudaError_t err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            printf("CUDA test kernel error: %s\n", cudaGetErrorString(err));
        }

        // Allocate host memory for matrices
        // Allow for matrices to come from device -> host
        // Free up device memory
        int *hostScoringMatrix = new int[(referenceLength+1) * (queryLength+1)];
        direction *hostBacktrackMatrix = new direction[(referenceLength+1) * (queryLength+1)];

        cudaMemcpy(hostScoringMatrix, deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(hostBacktrackMatrix, deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction), cudaMemcpyDeviceToHost);

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