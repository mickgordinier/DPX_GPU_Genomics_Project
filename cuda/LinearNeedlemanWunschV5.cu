#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"
#include "../c++/timing.h"
#include "pthread.h"

// Blocks are 1D with a size of the 32 threads (For 1 warp)
#define BLOCK_SIZE 32
#define BATCH_SIZE 1000

//#define BACKTRACK_BATCH_SIZE 10
//#define BACKTRACK_BATCH_SIZE 20
//#define BACKTRACK_BATCH_SIZE 100
#define BACKTRACK_BATCH_SIZE 200

//#define DEBUG
// Defining this will test all of the sequences in the input file
#define TEST_ALL

/*
    THINGS TO CONSIDER FOR OPTIMIZATION
    1. Complete removal of the scoring matrix altogether (Use of warp shuffling and shared memory)
    2. Using 16x2 DPX instructions to have a thread work on 2 cells concurrently

*/

struct thread_arg{
    int threadID;
    int pairNum;
    directionMain ** backtrackingMemos;
    char* sequences;
    seqPair* allSequenceInfo;
    int *scores;
};

void* threadBacktrackLNW(void *tmp){
	struct thread_arg* threadData = (struct thread_arg*)tmp;
    int threadID                = threadData->threadID;
    int pairNum                 = threadData->pairNum;
    directionMain ** backtrackingMemos       = threadData->backtrackingMemos;
    char* sequences             = threadData->sequences;
    seqPair* allSequenceInfo    = threadData->allSequenceInfo;
    int* scores                 = threadData->scores;
    
    #ifdef DEBUG
    printLock();
    printf("Thread ID: %d\n", threadID);
    fflush(stdout);
    printUnlock();
    #endif
    
    for(int i = 0; i < BACKTRACK_BATCH_SIZE; i++){
    
        const char *reference = &sequences[allSequenceInfo[pairNum + i].referenceIdx];
        const char *query = &sequences[allSequenceInfo[pairNum + i].queryIdx];
        
        const int referenceLength = allSequenceInfo[pairNum + i].referenceSize;
        const int queryLength = allSequenceInfo[pairNum + i].querySize;

        // Backtrrack
        backtrackMultiNW(backtrackingMemos[i], reference, referenceLength, query, queryLength, pairNum + i, scores[i]);
        
        // Free a matrix
        free(backtrackingMemos[i]);
    }
    // Free up memory
    free(backtrackingMemos); // Free entire matrix
    free(scores); // Free all scores
}

// NEEDLEMAN WUNSCH BASELINE KERNEL
// SUPPORTS KERNELS WITH THREADS LESS THAN QUERY LENGTH
// KERNEL IS LAUNCHED WITH MULTIPLE BLOCKS
// EACH BLOCK HANDLES A SEQUENCE PAIR

__global__ void 
needleman_wunsch_kernel(
    int **batchScoringMatrices, directionMain **batchBacktrackMatrices,
    const char *allSequences, const seqPair *allSequenceInfo,
    const int matchWeight, const int mismatchWeight, const int gapWeight,
    const int startingSequenceIdx)
{

    const int threadCount = blockDim.x;
    const int tid = threadIdx.x;

    // We are launching multiple blocks, each of a warp of threads
    // Each block handles their own sequence alignment
    // We index into the array to obtain the strings and length
    
    int *scoringMatrix = batchScoringMatrices[blockIdx.x];
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

        int **tempDeviceScoringMatrices = (int **)malloc(BATCH_SIZE * sizeof(int*));
        directionMain **tempDeviceBacktrackMatrices = (directionMain **)malloc(BATCH_SIZE * sizeof(directionMain*));

        int **deviceScoringMatrices;
        directionMain **deviceBacktrackMatrices;

        handleErrs(
            cudaMalloc(&deviceScoringMatrices, BATCH_SIZE * sizeof(int*)),
            "FAILED TO ALLOCATE MEMORY TO deviceScoringMatrices\n"
        );

        handleErrs(
            cudaMalloc(&deviceBacktrackMatrices, BATCH_SIZE * sizeof(directionMain*)),
            "FAILED TO ALLOCATE MEMORY TO deviceBacktrackMatrices\n"
        );
        memalloc_time += get_time() - start_memalloc;

        // Run the kernel on every sequence
        for(size_t sequenceIdx = 0; sequenceIdx < fileInfo.numPairs; sequenceIdx+=BATCH_SIZE){
            
            start_memalloc = get_time();
            for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; ++i) {
                
                const int queryLength = allSequenceInfo[i].querySize;
                const int referenceLength = allSequenceInfo[i].referenceSize;
                
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

                tempDeviceScoringMatrices[i-sequenceIdx] = deviceScoringMatrix;
                tempDeviceBacktrackMatrices[i-sequenceIdx] = deviceBacktrackMatrix;
            }

            handleErrs(
                cudaMemcpy(deviceScoringMatrices, tempDeviceScoringMatrices, BATCH_SIZE * sizeof(int*), cudaMemcpyHostToDevice),
                "FAILED TO COPY MEMORY FOR deviceScoringMatrices\n"
            );

            handleErrs(
                cudaMemcpy(deviceBacktrackMatrices, tempDeviceBacktrackMatrices, BATCH_SIZE * sizeof(directionMain*), cudaMemcpyHostToDevice),
                "FAILED TO COPY MEMORY FOR deviceBacktrackMatrices\n"
            );
            memalloc_time += get_time() - start_memalloc;

            uint64_t start_kernel = get_time();
            // Need to launch kernel
            needleman_wunsch_kernel<<<BATCH_SIZE, BLOCK_SIZE>>>(
                deviceScoringMatrices, deviceBacktrackMatrices,
                deviceSequences, deviceAllSequenceInfo,
                matchWeight, mismatchWeight, gapWeight,
                sequenceIdx
            );
            
            // Wait for kernel to finish
            handleErrs(
                cudaDeviceSynchronize(),
                "SYNCHRONIZATION FAILED\n"
            );
            kernel_time += get_time() - start_kernel;

            // Create all of our backtracking threads
            start_memalloc = get_time();
            pthread_t           *threads    = (pthread_t*)malloc(sizeof(pthread_t)*BATCH_SIZE/BACKTRACK_BATCH_SIZE);
            struct thread_arg   *threadArgs = (thread_arg*)malloc(sizeof(thread_arg)*BATCH_SIZE/BACKTRACK_BATCH_SIZE);
            memalloc_time += get_time() - start_memalloc;
            for (int pairBaseIdx = sequenceIdx; pairBaseIdx < sequenceIdx+BATCH_SIZE; pairBaseIdx += BACKTRACK_BATCH_SIZE) {
                // Init thread arguments
                int threadID =  (pairBaseIdx - sequenceIdx) / BACKTRACK_BATCH_SIZE;
                struct thread_arg threadArg;
                threadArg.threadID        = threadID;
                threadArg.pairNum         = pairBaseIdx;
                threadArg.sequences       = sequences;
                threadArg.allSequenceInfo = allSequenceInfo;
                
                start_memalloc = get_time();
                // Create array of scores for a single thread
                threadArg.scores = (int*)malloc(sizeof(int) * BACKTRACK_BATCH_SIZE);

                // Create array of backtracking matrices for a single thread
                threadArg.backtrackingMemos = (directionMain**)malloc(sizeof(directionMain*) * BACKTRACK_BATCH_SIZE);
                
                // Copy a subset of the scores and backtracking matrices back from cuda
                for(int pairIdx = pairBaseIdx; pairIdx < pairBaseIdx + BACKTRACK_BATCH_SIZE; pairIdx++){
                    const char *queryString = &sequences[allSequenceInfo[pairIdx].queryIdx];
                    const char *referenceString = &sequences[allSequenceInfo[pairIdx].referenceIdx];
                    
                    const int queryLength = allSequenceInfo[pairIdx].querySize;
                    const int referenceLength = allSequenceInfo[pairIdx].referenceSize;

                    // Get max score
                    int *hostScoringMatrix = (int*)malloc(sizeof(int) * (referenceLength+1) * (queryLength+1));
                    handleErrs(
                        cudaMemcpy(hostScoringMatrix, tempDeviceScoringMatrices[pairIdx-sequenceIdx], (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost),
                        "FAILED TO COPY SCORING MATRIX FROM DEVICE --> HOST\n"
                    );
                    threadArg.scores[pairIdx - pairBaseIdx] = hostScoringMatrix[(referenceLength + 1) * (queryLength + 1) - 1];
                    cudaFree(tempDeviceScoringMatrices[pairIdx-sequenceIdx]);
                    free(hostScoringMatrix);

                    // Get backtack matrix
                    threadArg.backtrackingMemos[pairIdx - pairBaseIdx] = (directionMain*)malloc(sizeof(directionMain)*(referenceLength+1) * (queryLength+1));
                    handleErrs(
                        cudaMemcpy(threadArg.backtrackingMemos[pairIdx - pairBaseIdx], tempDeviceBacktrackMatrices[pairIdx-sequenceIdx], (referenceLength+1) * (queryLength+1) * sizeof(directionMain), cudaMemcpyDeviceToHost),
                        "FAILED TO COPY BACKTRACK MATRIX FROM DEVICE --> HOST\n"
                    );
                    cudaFree(tempDeviceBacktrackMatrices[pairIdx-sequenceIdx]);
                    
                }
                memalloc_time += get_time() - start_memalloc;
                threadArgs[threadID] = threadArg;
                
                // Backtrack create a thread
                uint64_t start_backtrack = get_time();
                #ifdef DEBUG
                    printf("Creating threadID: %d\n", threadID);
                #endif
                int ret = pthread_create(&threads[threadID], NULL, threadBacktrackLNW, (void*) &threadArgs[threadID]);
                if(ret != 0){
                    printf("Failed to create threadID: %d", threadID);
                    exit(1);
                }
                backtracking_time += get_time() - start_backtrack;
            }

            // Wait for each thread to finish
            uint64_t start_backtrack = get_time();
            for(int threadID = 0; threadID < BATCH_SIZE/BACKTRACK_BATCH_SIZE; threadID++){
                #ifdef DEBUG
                    printLock();
                    printf("Attempting to wait for thread thread: %d\n", threadID);
                    fflush(stdout);
                    printUnlock();
                #endif
                int ret;
                ret = pthread_join(threads[threadID], NULL);
                if(ret != 0){
                    printf("Failde to join to threadID: %d", threadID);
                    exit(1);
                }
            }
            delete[] threads;
            delete[] threadArgs;
            backtracking_time += get_time() - start_backtrack;
        }

        cudaFree(deviceSequences);
        cudaFree(deviceAllSequenceInfo);

        free(tempDeviceScoringMatrices);
        free(tempDeviceBacktrackMatrices);

        cudaFree(deviceScoringMatrices);
        cudaFree(deviceBacktrackMatrices);
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
        uint64_t start_backtrack = get_time();
        backtrackNW(hostBacktrackMatrix, referenceString, referenceLength, queryString, queryLength);
        backtracking_time += get_time() - start_backtrack;
        
        // Free data arrays
        delete[] hostScoringMatrix;
        delete[] hostBacktrackMatrix;
    #endif

    uint64_t elapsed_time = get_elapsed_time();
    printf("Elapsed time (usec): %lld\n", elapsed_time);
    printf("Elapsed time kernel (usec): %lld\n", kernel_time);
    printf("Elapsed time backtracking (usec): %lld\n", backtracking_time);
    printf("Elapsed time memallocing (usec): %lld\n", memalloc_time);
    printf("Elapsed time sum (usec): %lld\n",kernel_time + backtracking_time + memalloc_time);

    // Cleanup
    printf("Cleaning up\n");
    cleanupParsedFile(allSequenceInfo, sequences);
}