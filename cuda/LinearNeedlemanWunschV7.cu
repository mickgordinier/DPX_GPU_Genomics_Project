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
#define BACKTRACK_BATCH_SIZE 100

//#define DEBUG
// Defining this will test all of the sequences in the input file
#define TEST_ALL

int numMemCpyDirections = 0;
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
        backtrackMultiNW(backtrackingMemos[(threadID*BACKTRACK_BATCH_SIZE) + i], reference, referenceLength, query, queryLength, pairNum + i, scores[(threadID*BACKTRACK_BATCH_SIZE) + i]);
        
        // Free a matrix
        delete[] backtrackingMemos[(threadID*BACKTRACK_BATCH_SIZE) + i];
    }
}

// NEEDLEMAN WUNSCH BASELINE KERNEL
// SUPPORTS KERNELS WITH THREADS LESS THAN QUERY LENGTH
// KERNEL IS LAUNCHED WITH MULTIPLE BLOCKS
// EACH BLOCK HANDLES A SEQUENCE PAIR

__global__ void 
needleman_wunsch_kernel(
    int *similarityScores,
    directionMain **batchBacktrackMatrices,
    const char *allSequences, const seqPair *allSequenceInfo,
    const int matchWeight, const int mismatchWeight, const int gapWeight,
    const int startingSequenceIdx)
{

    const int threadCount = blockDim.x;
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
        elementIdx += threadCount;
    }

    // Initialize the left col
    // NOT Writing in DRAM burst (slower)
    elementIdx = tid;
    while(elementIdx < numRows) {
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

    int leftDiag = gapWeight*tid;
    int left = gapWeight*(tid+1);
    int up = gapWeight*(tid+1); 
    int tmp_left = left;

    // Going through all of the rows each thread has to do
    for (int stripeStart = 1; stripeStart < numRows; stripeStart+=threadCount){

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

        
        // Allocate space for similarity scores and backtracking matrices
        int *hostSimilarityScores       = new int[BATCH_SIZE];
        int *deviceSimilarityScores;
        directionMain **tempDeviceBacktrackMatrices = new directionMain *[BATCH_SIZE];
        directionMain **deviceBacktrackMatrices;
        
        handleErrs(
            cudaMalloc(&deviceSimilarityScores, BATCH_SIZE * sizeof(int)),
            "FAILED TO ALLOCATE MEMORY TO deviceBacktrackMatrices\n"
        );
        
        handleErrs(
            cudaMalloc(&deviceBacktrackMatrices, BATCH_SIZE * sizeof(directionMain*)),
            "FAILED TO ALLOCATE MEMORY TO deviceBacktrackMatrices\n"
        );


        // Allocate space for threads and threadArgs
        pthread_t           *threads    = new pthread_t[BATCH_SIZE/BACKTRACK_BATCH_SIZE];
        struct thread_arg   *threadArgs = new thread_arg[BATCH_SIZE/BACKTRACK_BATCH_SIZE];
        directionMain      **hostBacktrackingMemos = new directionMain *[BATCH_SIZE];

        memalloc_time += get_time() - start_memalloc;
        
        // Main for loop to run the kernel on groups of [BATCH_SIZE] sequences
        for(size_t sequenceIdx = 0; sequenceIdx < fileInfo.numPairs; sequenceIdx+=BATCH_SIZE){
            
            start_memalloc = get_time();
            
            int largestReferenceLength = 0;

            // Allocate space for all of the backtracking matrices
            for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; ++i) {
                
                const int queryLength = allSequenceInfo[i].querySize;
                const int referenceLength = allSequenceInfo[i].referenceSize;
                
                largestReferenceLength = max(largestReferenceLength, referenceLength);
                
                directionMain *deviceBacktrackMatrix;
                
                handleErrs(
                    cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(directionMain)),
                    "FAILED TO ALLOCATE MEMORY TO BACKTRACK MATRIX\n"
                );
                
                // Add new backtracking matrix to our array of matrices
                tempDeviceBacktrackMatrices[i-sequenceIdx] = deviceBacktrackMatrix;
            }
            
            // Copy over the pointers to our array of matrices
            handleErrs(
                cudaMemcpy(deviceBacktrackMatrices, tempDeviceBacktrackMatrices, BATCH_SIZE * sizeof(directionMain*), cudaMemcpyHostToDevice),
                "FAILED TO COPY MEMORY FOR deviceBacktrackMatrices\n"
            );
            memalloc_time += get_time() - start_memalloc;
            
            uint64_t start_kernel = get_time();
            // Need to launch kernel
            int smem_size = (largestReferenceLength + 1) * sizeof(int);
            needleman_wunsch_kernel<<<BATCH_SIZE, BLOCK_SIZE, smem_size>>>(
                deviceSimilarityScores,
                deviceBacktrackMatrices,
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
            
            // Get matrix scores and backtracking matrices back
            start_memalloc = get_time();
            handleErrs(
                cudaMemcpy(hostSimilarityScores, deviceSimilarityScores, BATCH_SIZE * sizeof(int), cudaMemcpyDeviceToHost),
                "FAILED TO COPY SIMILARITY SCORES FROM DEVICE --> HOST\n"
            );
            // // Copy backtracking matrices back from cuda into hostBacktrackingMemos
            // for(int pairIdx = sequenceIdx; pairIdx < sequenceIdx + BATCH_SIZE; pairIdx++){
            //     const char *queryString = &sequences[allSequenceInfo[pairIdx].queryIdx];
            //     const char *referenceString = &sequences[allSequenceInfo[pairIdx].referenceIdx];

            //     const int queryLength = allSequenceInfo[pairIdx].querySize;
            //     const int referenceLength = allSequenceInfo[pairIdx].referenceSize;

            //     // Get backtack matrix
            //     directionMain *backtrackMatrix = new directionMain[(referenceLength+1) * (queryLength+1)];
            //     handleErrs(
            //             cudaMemcpy(backtrackMatrix, tempDeviceBacktrackMatrices[pairIdx-sequenceIdx], (referenceLength+1) * (queryLength+1) * sizeof(directionMain), cudaMemcpyDeviceToHost),
            //         "FAILED TO COPY BACKTRACK MATRIX FROM DEVICE --> HOST\n"
            //     );
            //     cudaFree(tempDeviceBacktrackMatrices[pairIdx-sequenceIdx]);

            //     hostBacktrackingMemos[pairIdx - sequenceIdx] = backtrackMatrix;
            // }
            memalloc_time += get_time() - start_memalloc;
            
            // Make batch of threads
            for (int pairBaseIdx = sequenceIdx; pairBaseIdx < sequenceIdx+BATCH_SIZE; pairBaseIdx += BACKTRACK_BATCH_SIZE) {
                // Init thread arguments
                int threadID =  (pairBaseIdx - sequenceIdx) / BACKTRACK_BATCH_SIZE;
                struct thread_arg threadArg;
                threadArg.threadID        = threadID;
                threadArg.pairNum         = pairBaseIdx;
                threadArg.sequences       = sequences;
                threadArg.allSequenceInfo = allSequenceInfo;
                
                // Point the thread to the scoring array
                threadArg.scores = hostSimilarityScores;

                // Point thread to the backtracking matrices
                threadArg.backtrackingMemos = hostBacktrackingMemos;
                
                threadArgs[threadID] = threadArg;


                // // Copy backtracking matrices back from cuda into hostBacktrackingMemos
                //start_memalloc = get_time();
                for(int pairIdx = pairBaseIdx; pairIdx < pairBaseIdx + BACKTRACK_BATCH_SIZE; pairIdx++){
                    const char *queryString = &sequences[allSequenceInfo[pairIdx].queryIdx];
                    const char *referenceString = &sequences[allSequenceInfo[pairIdx].referenceIdx];

                    const int queryLength = allSequenceInfo[pairIdx].querySize;
                    const int referenceLength = allSequenceInfo[pairIdx].referenceSize;

                    // Get backtack matrix
                    start_memalloc = get_time();
                    directionMain *backtrackMatrix = new directionMain[(referenceLength+1) * (queryLength+1)];
                    handleErrs(
                            cudaMemcpy(backtrackMatrix, tempDeviceBacktrackMatrices[pairIdx-sequenceIdx], (referenceLength+1) * (queryLength+1) * sizeof(directionMain), cudaMemcpyDeviceToHost),
                        "FAILED TO COPY BACKTRACK MATRIX FROM DEVICE --> HOST\n"
                    );
                    cudaFree(tempDeviceBacktrackMatrices[pairIdx-sequenceIdx]);
                    memalloc_time += get_time() - start_memalloc;

                    hostBacktrackingMemos[pairIdx - sequenceIdx] = backtrackMatrix;
                    numMemCpyDirections += 1;
                }
                //memalloc_time += get_time() - start_memalloc;

                // Create backtracking threads
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
                    printf("Attempting to wait for threadID: %d\n", threadID);
                    fflush(stdout);
                    printUnlock();
                #endif
                int ret;
                ret = pthread_join(threads[threadID], NULL);
                if(ret != 0){
                    printf("Failed to join to threadID: %d", threadID);
                    exit(1);
                }
            }
            backtracking_time += get_time() - start_backtrack;
        }
        
        cudaFree(deviceSequences);
        cudaFree(deviceAllSequenceInfo);
        
        // free(tempDeviceScoringMatrices);
        delete[] tempDeviceBacktrackMatrices;
        delete[] hostSimilarityScores;
        delete[] hostBacktrackingMemos;
        
        // cudaFree(deviceScoringMatrices);
        cudaFree(deviceBacktrackMatrices);
        cudaFree(deviceSimilarityScores);

        delete[] threads;
        delete[] threadArgs;
    #endif
    uint64_t elapsed_time = get_elapsed_time();
    printf("Elapsed time (usec): %lld\n", elapsed_time);
    printf("Elapsed kernel time (usec): %lld\n", kernel_time);
    printf("Elapsed backtracking time (usec): %lld\n", backtracking_time);
    printf("Elapsed memallocing time (usec): %lld\n", memalloc_time);
    printf("Elapsed time sum (usec): %lld\n",kernel_time + backtracking_time + memalloc_time);
    printf("Num mem copys: %lld\n",numMemCpyDirections);

    // Cleanup
    printf("Cleaning up\n");
    cleanupParsedFile(allSequenceInfo, sequences);
}