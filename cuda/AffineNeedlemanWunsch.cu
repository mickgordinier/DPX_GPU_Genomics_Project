#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"
#include "../c++/timing.h"

// Blocks are 1D with a size of the 32 threads (For 1 warp)
#define BLOCK_SIZE 32

// Defining this will test all of the sequences in the input file
#define TEST_ALL

/*
    THINGS TO CONSIDER FOR OPTIMIZATION
    1. Complete removal of the scoring matrix altogether (Use of warp shuffling and shared memory)
    2. Using 16x2 DPX instructions to have a thread work on 2 cells concurrently

*/

// NEEDLEMAN WUNSCH BASELINE KERNEL

__global__ void 
affine_needleman_wunsch_kernel(
    int *scoringMatrix, directionMain *backtrackScoringMatrix,
    int *queryDeletionMatrix, directionIndel *backtrackDeletionMatrix,
    int *queryInsertionMatrix, directionIndel *backtrackInsertionMatrix,
    const char *queryString, const char *referenceString,
    const int queryLength, const int referenceLength,
    const int matchWeight, const int mismatchWeight, 
    const int gapOpenWeight, const int gapExtendWeight)
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
        scoringMatrix[elementIdx] = gapOpenWeight + (elementIdx * gapExtendWeight);
        backtrackScoringMatrix[elementIdx] = QUERY_INSERTION;
        elementIdx += threadCount;
    }

    // Initialize the left col
    // NOT Writing in DRAM burst (slower)
    elementIdx = tid;
    while(elementIdx < numRows) {
        scoringMatrix[elementIdx*numCols] = gapOpenWeight + (elementIdx * gapExtendWeight);
        backtrackScoringMatrix[elementIdx*numCols] = QUERY_DELETION;
        elementIdx += threadCount;
    }

    if (tid == 0) {
        backtrackScoringMatrix[0] = NONE_MAIN;
        backtrackDeletionMatrix[0] = NONE_INDEL;
        backtrackInsertionMatrix[0] = NONE_INDEL;
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

    // Each thread needs to iterate through the loop to be able to make the __syncthreads() call
    // All threads need to be able to reach the __syncthreads() call
    const int differentRows = ((numRows - 1) / BLOCK_SIZE) + 1;

    // Every thread gets a row and char
    int rowIdx = tid + 1;
    char queryChar = queryString[rowIdx - 1];

    int cell00Idx;
    int cell01Idx;
    int cell10Idx;
    int cell11Idx;

    for (int rowLoopIdx = 0; rowLoopIdx < differentRows; ++rowLoopIdx) {

        // If the thread in the warp is outside the matrix, wait for the other threads
        if (rowIdx < numRows) {

            // Each later thread must wait for the previous thread
            int adjCol = 1 - tid;

            // Setup cell indices once a thread can start executing
            cell00Idx = (rowIdx-1)*numCols;
            cell01Idx = (rowIdx-1)*numCols + 1;
            cell10Idx = rowIdx*numCols;
            cell11Idx = rowIdx*numCols + 1;
            
            // Each thread must go through the whole row
            // BUT, there is an adjustment that each thread must wait for
            for (int colIdx = 1; colIdx < (numCols+BLOCK_SIZE); ++colIdx) {

                // Main cell updating
                if((adjCol > 0) && (adjCol < numCols)){

                    bool pred;

                    // Handling scores of performing an query deletion at the end
                    // Calculating best score of either creating or extending the deletion gap
                    if (rowIdx == 1) {
                        // PROBABLY CAN HANDLE ROW 1 DURING INITIALIZATION PHASE
                        // Always assuming just opening new gap
                        queryDeletionMatrix[cell11Idx] = scoringMatrix[cell01Idx] + gapOpenWeight + gapExtendWeight;
                        backtrackDeletionMatrix[cell11Idx] = GAP_OPEN;
                    } else {
                        queryDeletionMatrix[cell11Idx] = __vibmax_s32(
                            scoringMatrix[cell01Idx] + gapOpenWeight + gapExtendWeight,  // Opening new gap at the end
                            queryDeletionMatrix[cell01Idx] + gapExtendWeight,            // Extending current gap at end
                            &pred
                        );
                        backtrackDeletionMatrix[cell11Idx] = pred ? GAP_OPEN : GAP_EXTEND;
                    }

                    // Handling scores of performing an query insertion at the end
                    // Calculating best score of either creating or extending the insertion gap
                    if (colIdx == 1) {
                        // PROBABLY CAN HANDLE COL 1 DURING INITIALIZATION PHASE
                        // Always assuming just opening new gap
                        queryInsertionMatrix[cell11Idx] = scoringMatrix[cell10Idx] + gapOpenWeight + gapExtendWeight;
                        backtrackInsertionMatrix[cell11Idx] = GAP_OPEN;
                    } else {
                        queryInsertionMatrix[cell11Idx] = __vibmax_s32(
                            scoringMatrix[cell10Idx] + gapOpenWeight + gapExtendWeight,  // Opening new gap at the end
                            queryInsertionMatrix[cell10Idx] + gapExtendWeight,           // Extending current gap at end
                            &pred
                        );
                        backtrackInsertionMatrix[cell11Idx] = pred ? GAP_OPEN : GAP_EXTEND;
                    }
                    
                    char referenceChar = referenceString[adjCol - 1];
                    directionMain cornerDirection;

                    // Determine if match
                    bool isMatch = (queryChar == referenceChar);
                    cornerDirection = isMatch ? MATCH : MISMATCH;
    
                    // Get all the possible scores
                    int matchMismatchScore = isMatch ? scoringMatrix[cell00Idx] + matchWeight : scoringMatrix[cell00Idx] + mismatchWeight;
            
                    int largestScore = __vibmax_s32(queryDeletionMatrix[cell11Idx], matchMismatchScore, &pred);
                    if (pred) cornerDirection = QUERY_DELETION;
                    
                    largestScore = __vibmax_s32(queryInsertionMatrix[cell11Idx], largestScore, &pred);
                    if (pred) cornerDirection = QUERY_INSERTION;

                    scoringMatrix[cell11Idx] = largestScore;
                    backtrackScoringMatrix[cell11Idx] = cornerDirection;
                    
                    cell00Idx += 1;
                    cell01Idx += 1;
                    cell10Idx += 1;
                    cell11Idx += 1;
                }

                ++adjCol;

            } // end

        } // end if

        // All previous threads must finish before moving onto the next row
        __syncthreads();

        rowIdx += BLOCK_SIZE;
        queryChar = queryString[rowIdx - 1];

    } // end for

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
    int gapOpenWeight   = -3;
    int gapExtendWeight = -1;
    if(strcmp(argv[1], "-pairs") == 0) {
        pairFileName = argv[2];
    }
    if(argc > 3 && strcmp(argv[3], "-match") == 0) {
        matchWeight = atoi(argv[4]);
    }
    if(argc > 5 && strcmp(argv[5], "-mismatch") == 0) {
        mismatchWeight = atoi(argv[6]);
    }
    if(argc > 7 && strcmp(argv[7], "-open") == 0) {
        gapOpenWeight = atoi(argv[8]);
    }
    if(argc > 9 && strcmp(argv[9], "-extend") == 0) {
        gapExtendWeight = atoi(argv[10]);
    }

    // Parse input file
    printf("Parsing input file: %s\n", pairFileName);
    inputInfo fileInfo;
    seqPair* sequenceIdxs;
    char* sequences;
    fileInfo = parseInput(pairFileName, sequenceIdxs, sequences);
    printf("Num Pairs: %d\n\n", fileInfo.numPairs);

    // Start timer
    uint64_t kernel_time = 0;
    uint64_t memalloc_time = 0;
    uint64_t backtracking_time = 0;
    uint64_t start_time = start_timer();
        
    // Copy over the sequences
    char* deviceSequences;
    uint64_t start_memalloc = get_time();
    handleErrs(
        cudaMalloc(&deviceSequences, (fileInfo.numBytes) * sizeof(char)),
        "FAILED TO ALLOCATE MEMORY FOR ALL SEQUENCES\n"
    );

    handleErrs(
        cudaMemcpy(deviceSequences, sequences, (fileInfo.numBytes) * sizeof(char), cudaMemcpyHostToDevice),
        "FAILED TO COPY MEMORY FOR ALL SEQUENCES\n"
    );

    memalloc_time += get_time() - start_memalloc;

    // Run the kernel on every sequence
    for(size_t i = 0; i < fileInfo.numPairs; i++){
        start_memalloc = get_time();

        char *referenceString = &sequences[sequenceIdxs[i].referenceIdx];
        char *queryString = &sequences[sequenceIdxs[i].queryIdx];

        int referenceLength = strlen(referenceString);
        int queryLength = strlen(queryString);

        int *deviceScoringMatrix;
        directionMain *deviceBacktrackMatrix;

        int *deviceQueryDeletionMatrix;
        directionIndel *deviceDeletionBacktrackMatrix;

        int *deviceQueryInsertionMatrix;
        directionIndel *deviceInsertionBacktrackMatrix;

        handleErrs(
            cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)),
            "FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n"
        );

        handleErrs(
            cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(directionMain)),
            "FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n"
        );

        handleErrs(
            cudaMalloc(&deviceQueryDeletionMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)),
            "FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n"
        );

        handleErrs(
            cudaMalloc(&deviceDeletionBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(directionIndel)),
            "FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n"
        );

        handleErrs(
            cudaMalloc(&deviceQueryInsertionMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int)),
            "FAILED TO ALLOCATE MEMORY TO SCORING MATRIX\n"
        );

        handleErrs(
            cudaMalloc(&deviceInsertionBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(directionIndel)),
            "FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n"
        );

        memalloc_time += get_time() - start_memalloc;
        uint64_t start_kernel = get_time();

        // Need to launch kernel
        affine_needleman_wunsch_kernel<<<1, BLOCK_SIZE>>>(
            deviceScoringMatrix, deviceBacktrackMatrix,
            deviceQueryDeletionMatrix, deviceDeletionBacktrackMatrix,
            deviceQueryInsertionMatrix, deviceInsertionBacktrackMatrix,
            deviceSequences + sequenceIdxs[i].queryIdx, deviceSequences + sequenceIdxs[i].referenceIdx, 
            sequenceIdxs[i].querySize, sequenceIdxs[i].referenceSize, 
            matchWeight, mismatchWeight, gapOpenWeight, gapExtendWeight
        );
        
        // Wait for kernel to finish
        handleErrs(
            cudaDeviceSynchronize(),
            "SYNCHRONIZATION FAILED\n"
        );

        kernel_time += get_time() - start_kernel;
        start_memalloc = get_time();

        // Copy the matrices back over
        int *hostScoringMatrix = new int[(referenceLength+1) * (queryLength+1)];
        directionMain *hostBacktrackMatrix = new directionMain[(referenceLength+1) * (queryLength+1)];

        int *hostQueryDeletionMatrix = new int[(referenceLength+1) * (queryLength+1)];
        directionIndel *hostDeletionBacktrackMatrix = new directionIndel[(referenceLength+1) * (queryLength+1)];

        int *hostQueryInsertionMatrix = new int[(referenceLength+1) * (queryLength+1)];
        directionIndel *hostInsertionBacktrackMatrix = new directionIndel[(referenceLength+1) * (queryLength+1)];

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

        cudaFree(deviceQueryDeletionMatrix);
        cudaFree(deviceDeletionBacktrackMatrix);

        cudaFree(deviceQueryInsertionMatrix);
        cudaFree(deviceInsertionBacktrackMatrix);

        memalloc_time += get_time() - start_memalloc;

        uint64_t start_backtrack = get_time();
        // Backtrack matrices
        printf("%d | %d\n", i, hostScoringMatrix[(referenceLength + 1) * (queryLength + 1) - 1]);
        backtrackANW(hostBacktrackMatrix, hostInsertionBacktrackMatrix, hostDeletionBacktrackMatrix, referenceString, referenceLength, queryString, queryLength);
        backtracking_time += get_time() - start_backtrack;

        // Free data arrays
        delete[] hostScoringMatrix;
        delete[] hostBacktrackMatrix;

        delete[] deviceQueryDeletionMatrix;
        delete[] deviceDeletionBacktrackMatrix;

        delete[] deviceQueryInsertionMatrix;
        delete[] deviceInsertionBacktrackMatrix;
    }

    cudaFree(deviceSequences);

    uint64_t elapsed_time = get_elapsed_time();
    printf("Elapsed time (usec): %lld\n", elapsed_time);
    printf("Elapsed kernel time (usec): %lld\n", kernel_time);
    printf("Elapsed backtracking time (usec): %lld\n", backtracking_time);
    printf("Elapsed memallocing time (usec): %lld\n", memalloc_time);
    printf("Elapsed time sum (usec): %lld\n",kernel_time + backtracking_time + memalloc_time);

    // Cleanup
    printf("Cleaning up\n");
    cleanupParsedFile(sequenceIdxs, sequences);
}