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

/* TODO: add band lmao */

/* 
NOTES
-   each block handles a query-reference pair
-   each thread operates on a single element in a diagonal
-   threads iterate over the anti-diagonals of the matrix
*/
__global__ void banded_smith_waterman(int *scoringMatrix, direction *backtrackMatrix,
                        const char *queryString, const char *referenceString,
                        const int queryLength, const int referenceLength,
                        const int matchWeight, const int mismatchWeight, 
                        const int gapWeight,   const int threadCount){
    int tid = threadIdx.x;
    int block = blockIdx.x;
    int maxLen = max(queryLength, referenceLength);
    int minLen = min(queryLength, referenceLength);

    /* each block handles one query-reference pair */
    extern __shared__ char shared_mem[];
    char* shared_query = shared_mem;
    char* shared_ref = shared_query + queryLength;

    /* store scores for current and previous 2 diagonals */
    int * shared_diag_base = (int *)(shared_ref + referenceLength);
    int * current_diag = shared_diag_base;
    int * prev_diag = shared_diag_base + maxLen + 1;
    int * prev_prev_diag = shared_diag_base + 2*(maxLen + 1);

    __shared__ int max_row, max_col, max_score;

    /* load query & reference strings into shared memory */
    if (tid == 0){
        prev_diag[0] = 0;
        prev_prev_diag[0] = 0;
        current_diag[0] = 0;
        max_score = -1;
    }

    if (tid < queryLength){
        shared_query[tid] = queryString[tid];
    }

    if (tid < referenceLength){
        shared_ref[tid] = referenceString[tid];
    }

    __syncthreads();

    /* iterate over diagonals */
    for (int d = 2; d < minLen * 2 + 1; d++){
        
    }


}

int main(int argc, char *argv[]) {
    /* log Cuda details */
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

    /* ensure proper usage s*/
    if (argc < 2) {
		fprintf(stderr, "usage: main -pairs <InSeqFile> -match <matchWeight> -mismatch <mismatchWeight> -gap <gapWeight> \n");
		exit(EXIT_FAILURE);
    }
	
    /* get args */
    char *pairFileName;
    int matchWeight     = 3;
    int mismatchWeight  = -1;
    int gapWeight       = -2;
    int threadCount     = 32;
    // int band            = 16;
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
    // if (argc > 11 && strcmp(argv[11], "-band") == 0){
    //     band = atoi(argv[12]);
    // }

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