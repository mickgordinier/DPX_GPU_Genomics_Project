#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"

// Blocks are 1D with a size of the 32 threads (For 1 warp)
#define BLOCK_SIZE 32

// Defing this will test all of the sequences in the input file
#define TEST_ALL

/*
    THINGS TO CONSIDER FOR OPTIMIZATION
    1. Complete removal of the scoring matrix altogether (Use of warp shuffling and shared memory)
    2. Using 16x2 DPX instructions to have a thread work on 2 cells concurrently

*/

// NEEDLEMAN WUNSCH BASELINE KERNEL

__global__ void 
needleman_wunsch_kernel(
    int *scoringMatrix, directionMain *backtrackMatrix,
    const char *queryString, const char *referenceString,
    const int queryLength, const int referenceLength,
    const int matchWeight, const int mismatchWeight, 
    const int gapWeight)
{
    // We are only launching 1 block
    // Thus, each thread will only have a unique threadID that differentiates the threads
    const int tid = threadIdx.x;
    const int threadCount = blockDim.x;

    // The matrices are of size (queryLength + 1) * (referenceLength + 1)
    const int numRows = queryLength + 1;
    const int numCols = referenceLength + 1;

    /* --- (BEGIN) INITIALIZING THE SCORING MATRIX --- */

    // Used for when a thread has to iterate over more than one col/row
    int elementIdx;

    // Initialize the top row
    // Writing in DRAM burst for faster updating
    elementIdx = tid;
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

    // Each thread has to do ((numRows-1) / BLOCK_SIZE) different rows
    // During a row pass, each thread will be doing numCols amount of work
    // There is an additional numRows lag amount of waiting for staggering the threads
    const int numIterations = (numCols * ((numRows / BLOCK_SIZE) + 1)) + numRows;

    // Every thread gets a row and char
    int rowIdx = tid + 1;
    char queryChar = queryString[rowIdx - 1];

    // Each later thread must wait for the previous thread
    int adjCol = 1 - tid;

    // Setup cell indices once a thread can start executing
    int cell00Idx = (rowIdx-1)*numCols;
    int cell01Idx = (rowIdx-1)*numCols + 1;
    int cell10Idx = rowIdx*numCols;
    int cell11Idx = rowIdx*numCols + 1; 

    for (int iter = 0; iter < numIterations; ++iter) {
        
        if ((rowIdx < numRows) && (adjCol > 0)) {
            
            char referenceChar = referenceString[adjCol - 1];
            directionMain cornerDirection = NONE_MAIN;
            bool pred;
            
            // Determine if match
            bool isMatch = (queryChar == referenceChar);
            cornerDirection = isMatch ? MATCH : MISMATCH;

            // Get all the possible scores
            int matchMismatchScore = isMatch ? scoringMatrix[cell00Idx] + matchWeight : scoringMatrix[cell00Idx] + mismatchWeight;
            int queryDeletionScore = scoringMatrix[cell01Idx] + gapWeight;
            int queryInsertionScore = scoringMatrix[cell10Idx] + gapWeight;

            // Find the largest of the 3 scores
            // Utilizing DPX instructions for updating
            // pred = (queryDeletionScore >= matchMismatchScore)
            int largestScore;
            largestScore = __vibmax_s32(queryDeletionScore, matchMismatchScore, &pred);
            if (pred) cornerDirection = QUERY_DELETION;
            
            largestScore = __vibmax_s32(queryInsertionScore, largestScore, &pred);
            if (pred) cornerDirection = QUERY_INSERTION;

            // Update scoring matrix and incrementing pointers
            scoringMatrix[cell11Idx] = largestScore;
            backtrackMatrix[cell11Idx] = cornerDirection;
            cell00Idx += 1;
            cell01Idx += 1;
            cell10Idx += 1;
            cell11Idx += 1;

        }
        
        if (adjCol == numCols-1) {
            rowIdx += 32;
            queryChar = queryString[rowIdx - 1];
            adjCol = min(0, numCols-33);
            
            cell00Idx = (rowIdx-1)*numCols;
            cell01Idx = (rowIdx-1)*numCols + 1;
            cell10Idx = rowIdx*numCols;
            cell11Idx = rowIdx*numCols + 1; 
        }

        ++adjCol;
    }

    __syncthreads();

    /* --- (END) POPULATING THE SCORING MATRIX -- */
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
    // int threadCount     = 32;
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
    // if(argc > 9 && strcmp(argv[9], "-threads-per-alignment") == 0) {
    //     threadCount = atoi(argv[10]);
    // }

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
            directionMain *deviceBacktrackMatrix;

            handleErrs(
                cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)),
                "FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n"
            );
    
            handleErrs(
                cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(directionMain)),
                "FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n"
            );

            // Need to launch kernel
            needleman_wunsch_kernel<<<1, BLOCK_SIZE>>>(
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
            directionMain *hostBacktrackMatrix = new directionMain[(referenceLength+1) * (queryLength+1)];

            // Copy information back from device --> host
            handleErrs(
                cudaMemcpy(hostScoringMatrix, deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost),
                "FAILED TO COPY SCORING MATRIX FROM DEVICE --> HOST\n"
            );
            
            handleErrs(
                cudaMemcpy(hostBacktrackMatrix, deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(directionMain), cudaMemcpyDeviceToHost),
                "FAILED TO COPY BACKTRACK MATRIX FROM DEVICE --> HOST\n"
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
        // Launching a kernel with 1 block with threadCount threads to populate scoring matrix
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
        backtrackNW(hostBacktrackMatrix, referenceString, referenceLength, queryString, queryLength);
        
        // Free data arrays
        delete[] hostScoringMatrix;
        delete[] hostBacktrackMatrix;
    #endif

    // Cleanup
    printf("Cleaning up\n");
    cleanupParsedFile(sequenceIdxs, sequences);
}