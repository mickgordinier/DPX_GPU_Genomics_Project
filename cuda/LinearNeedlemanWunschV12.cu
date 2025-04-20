#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"
#include "../c++/timing.h"

// FOR PROFILING
#include <cupti_version.h>
#include <cupti.h>

// Blocks are 1D with a size of the 32 threads (For 1 warp)
#define BLOCK_SIZE 32
#define BATCH_SIZE 500

// Defining this will test all of the sequences in the input file
#define TEST_ALL

// CUPTI callback functions
void CUPTIAPI bufferRequested(uint8_t **buffer, size_t *size, size_t *maxNumRecords) {
    *size = 16 * 1024;
    *maxNumRecords = 0;
    *buffer = (uint8_t *)malloc(*size);
    if (*buffer == NULL) {
        printf("Error: out of memory\n");
        exit(1);
    }
}

void CUPTIAPI bufferCompleted(CUcontext ctx, uint32_t streamId, uint8_t *buffer, size_t size, size_t validSize) {
    CUpti_Activity *record = NULL;
    if (validSize > 0) {
        do {
            CUptiResult status = cuptiActivityGetNextRecord(buffer, validSize, &record);
            if (status == CUPTI_SUCCESS) {
                if (record->kind == CUPTI_ACTIVITY_KIND_KERNEL) {
                    CUpti_ActivityKernel4 *kernelRecord = (CUpti_ActivityKernel4 *)record;
                    printf("Kernel %s executed; Grid (%d, %d, %d), Block(%d, %d, %d), time: %lld\n",
                           kernelRecord->name,
                           kernelRecord->gridX, kernelRecord->gridY, kernelRecord->gridZ,
                           kernelRecord->blockX, kernelRecord->blockY, kernelRecord->blockZ,
                           kernelRecord->end - kernelRecord->start);
                }
                else if (record->kind == CUPTI_ACTIVITY_KIND_MEMSET) {
                    CUpti_ActivityMemset4  *memsetRecord  = (CUpti_ActivityMemset4 *)record;
                    printf("Memset: %llu bytes to device %llu in %llu ns\n",
                           memsetRecord->bytes, memsetRecord->deviceId,
                           memsetRecord->end - memsetRecord->start);
                }
                else if (record->kind == CUPTI_ACTIVITY_KIND_MEMCPY) {
                    CUpti_ActivityMemcpy2  *memcpyRecord  = (CUpti_ActivityMemcpy2 *)record;
                    printf("Memcpy: %llu bytes to device %llu in %llu ns\n",
                           memcpyRecord->bytes, memcpyRecord->deviceId,
                           memcpyRecord->end - memcpyRecord->start);
                }
                // else if (record->kind == CUPTI_ACTIVITY_KIND_GLOBAL_ACCESS) {
                //     // This is illustrative; replace with the actual structure type and data extraction
                //     CUpti_ActivityGlobalAccess3 *globalAccessRecord = (CUpti_ActivityGlobalAccess3 *)record;
                //     printf("Global access: executed %d\n", globalAccessRecord->executed);
                //     // Extract further details based on the actual structure
                // }
            }
            else if (status == CUPTI_ERROR_MAX_LIMIT_REACHED) break;
            else {
                printf("Error fetching outstanding records\n");
                break;
            }
        } while (1);
    }
    free(buffer);
}

/*
    THINGS TO CONSIDER FOR OPTIMIZATION
    1. Complete removal of the scoring matrix altogether (Use of warp shuffling and shared memory)
    2. Using 16x2 DPX instructions to have a thread work on 2 cells concurrently

*/

// NEEDLEMAN WUNSCH BASELINE KERNEL
// SUPPORTS KERNELS WITH THREADS LESS THAN QUERY LENGTH
// KERNEL IS LAUNCHED WITH MULTIPLE BLOCKS
// EACH BLOCK HANDLES A SEQUENCE PAIR

__global__ void 
needleman_wunsch_kernel(
    int *similarityScores,
    int *stringSpacing,
    directionMain *batchBacktrackMatrices,
    int *batchIndices,
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

    CUptiResult cuptiErr;
    CUpti_SubscriberHandle subscriber;

    // Setup CUPTI callback subscriber
    cuptiErr = cuptiSubscribe(&subscriber, (CUpti_CallbackFunc)bufferCompleted, NULL);
    if (cuptiErr != CUPTI_SUCCESS) {
        printf("Failed to subscribe CUPTI\n");
        exit(1);
    }

    cuptiErr = cuptiActivityRegisterCallbacks(bufferRequested, bufferCompleted);
    if (cuptiErr != CUPTI_SUCCESS) {
        printf("Failed to register CUPTI callbacks\n");
        exit(1);
    }

    // Enable activity types (e.g., CUPTI_ACTIVITY_KIND_KERNEL)
    cuptiErr = cuptiActivityEnable(CUPTI_ACTIVITY_KIND_KERNEL);
    if (cuptiErr != CUPTI_SUCCESS) {
        printf("Failed to enable CUPTI kernel activity\n");
        exit(1);
    }

    cuptiErr = cuptiActivityEnable(CUPTI_ACTIVITY_KIND_MEMCPY);
    if (cuptiErr != CUPTI_SUCCESS) {
        printf("Failed to enable CUPTI memcpy activity\n");
        exit(1);
    }

    cuptiErr = cuptiActivityEnable(CUPTI_ACTIVITY_KIND_MEMSET);
    if (cuptiErr != CUPTI_SUCCESS) {
        printf("Failed to enable CUPTI memset activity\n");
        exit(1);
    }

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

             /* memset to null bytes */
             handleErrs(
                cudaMemset(deviceBacktrackStringRet, 0, (stringLengthMax * 3) * BATCH_SIZE * sizeof(char)),
                "FAILED TO memset deviceMatricesAll\n"
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

                // ATTEMPT #1: BINARY SEARCH
                // int begin = ((stringLengthMax * 3) * (i - sequenceIdx));
                // int end = begin + stringLengthMax - 1;

                // int result = -1; 

                // while (begin <= end) {
                //     int middle = begin + (end - begin) / 2;

                //     if (hostBacktrackingStringRet[middle] == '\0') {
                //         begin = middle + 1;
                //     } else {
                //         result = middle;
                //         end = middle - 1;
                //     }
                // }

                // ATTEMPT #2: STUPID BINARY SEARCH
                // int spacing = ((stringLengthMax * 3) * (i - sequenceIdx));
                // int spacingIter = stringLengthMax / 4 - 1;
                
                // while (spacingIter > 0) {
                //     while (hostBacktrackingStringRet[spacing] == '\0' && (spacing < ((stringLengthMax * 3) * (i - sequenceIdx)) + stringLengthMax)) {
                //         spacing += spacingIter;
                //     }
                //     spacing -= spacingIter;

                //     spacingIter = spacingIter >> 1;
                // }
                // spacing++;

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
    #endif

    // Before disabling CUPTI activities or application exit
    cuptiErr = cuptiActivityFlushAll(0);
    if (cuptiErr != CUPTI_SUCCESS) {
        printf("Failed to flush CUPTI activities before exit\n");
    }

    // Clean up CUPTI
    cuptiErr = cuptiActivityDisable(CUPTI_ACTIVITY_KIND_KERNEL);
    if (cuptiErr != CUPTI_SUCCESS) {
        printf("Failed to disable CUPTI activity\n");
        exit(1);
    }

    cuptiErr = cuptiUnsubscribe(subscriber);
    if (cuptiErr != CUPTI_SUCCESS) {
        printf("Failed to unsubscribe CUPTI\n");
        exit(1);
    }


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