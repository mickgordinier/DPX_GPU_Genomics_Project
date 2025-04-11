#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"

// Blocks are 1D with a size of the maximum 1024 threads
#define BLOCK_SIZE 1024

/*
    THINGS TO CONSIDER FOR OPTIMIZATION
    1. Complete removal of the scoring matrix altogether (Use of warp shuffling and shared memory)
    2. Using 16x2 DPX instructions to have a thread work on 2 cells concurrently

*/

__global__ void test_kernel(int *scoringMatrix, direction *backtrackMatrix,
                        const char *queryString, const char *referenceString,
                        const int queryLength, const int referenceLength,
                        const int matchWeight, const int mismatchWeight, const int gapWeight){
    int tid = threadIdx.x;

    int numRows = queryLength + 1;
    int numCols = referenceLength + 1;
    // Initialize the matrix
    if (tid < numRows*numCols) {
        scoringMatrix[tid] = tid;
    }
    __syncthreads();

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
}

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


int main(int argc, char *argv[]) {

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
	
    // Get filename from args
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
    size_t numPairs;
    seqPair* sequenceIdxs;
    char* sequences;
    numPairs = parseInput(pairFileName, sequenceIdxs, sequences);

    // char *referenceString = &sequences[sequenceIdxs[0].referenceIdx];
    // char *queryString = &sequences[sequenceIdxs[0].queryIdx];
    char *referenceString = "GTCATGCAATAACG";
    char *queryString = "ATGCAATA";

    int referenceLength = strlen(referenceString);
    int queryLength = strlen(queryString);

    printf("Reference String: %s (Length: %d)\n", referenceString, referenceLength);
    printf("Query String: %s (Length: %d)\n", queryString, queryLength);
    printf("(MATCH WEIGHT, MISMATCH WEIGHT, GAP WEIGHT): (%d, %d, %d)\n\n", matchWeight, mismatchWeight, gapWeight);

    // Allocate device memory for matrices
    printf("[Allocating CUDA Memory]\n");
    int *deviceScoringMatrix;
    direction *deviceBacktrackMatrix;

    if (cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)) != cudaSuccess) {
        printf("FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n");
        exit(1);
    }

    if (cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction)) != cudaSuccess) {
        printf("FAILED TO ALLOCATE MEMORY TO BACKTRAK MATRIX\n");
        exit(1);
    }

    int *testDeviceScoringMatrix;
    direction *testDeviceBacktrackMatrix;
    if (cudaMalloc(&testDeviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)) != cudaSuccess) {
        printf("FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n");
        exit(1);
    }

    if (cudaMalloc(&testDeviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction)) != cudaSuccess) {
        printf("FAILED TO ALLOCATE MEMORY TO BACKTRAK MATRIX\n");
        exit(1);
    }
    

    // Need to launch kernel
    // Launching a kernel with 1 block with BLOCK_SIZE threads to populate scoring matrix
    needleman_wunsch_forward_pass_kernel<<<1, BLOCK_SIZE>>>(
        deviceScoringMatrix, deviceBacktrackMatrix,
        queryString, referenceString, queryLength, referenceLength, 
        matchWeight, mismatchWeight, gapWeight);
    test_kernel<<<1, BLOCK_SIZE>>>(
        testDeviceScoringMatrix, testDeviceBacktrackMatrix,
        queryString, referenceString, queryLength, referenceLength, 
        matchWeight, mismatchWeight, gapWeight);

    // Allocate host memory for matrices
    // Allow for matrices to come from device -> host
    // Free up device memory
    int *hostScoringMatrix = new int[(referenceLength+1) * (queryLength+1)];
    direction *hostBacktrackMatrix = new direction[(referenceLength+1) * (queryLength+1)];

    // Wait for kernel to finish
    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        printf("CUDA kernel error: %s\n", cudaGetErrorString(err));
    }

    cudaMemcpy(hostScoringMatrix, deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hostBacktrackMatrix, deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction), cudaMemcpyDeviceToHost);

    int *testHostScoringMatrix = new int[(referenceLength+1) * (queryLength+1)];
    direction *testHostBacktrackMatrix = new direction[(referenceLength+1) * (queryLength+1)];
    cudaMemcpy(testHostScoringMatrix, testDeviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(testHostBacktrackMatrix, testDeviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction), cudaMemcpyDeviceToHost);

    cudaFree(deviceScoringMatrix);
    cudaFree(deviceBacktrackMatrix);

    cudaFree(testDeviceScoringMatrix);
    cudaFree(testDeviceBacktrackMatrix);

    // Print Matrix
    printf("Test Scored Matrix\n");
    printMatrix(testHostScoringMatrix, referenceLength + 1, queryLength + 1);
    printf("Scored Matrix\n");
    printMatrix(hostScoringMatrix, referenceLength + 1, queryLength + 1);

    // Perform backtracking on host and print results
    printf("0 | %d", hostScoringMatrix[(referenceLength + 1) * (queryLength + 1) - 1]);
    backtrackNW(hostBacktrackMatrix, referenceString, referenceLength, queryString, queryLength);

    // Free data arrays
    printf("Cleaning up\n");
    cleanupParsedFile(sequenceIdxs, sequences);
    delete[] hostScoringMatrix;
    delete[] hostBacktrackMatrix;
}