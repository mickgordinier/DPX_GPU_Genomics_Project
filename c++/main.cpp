#include <cstring>
#include <cassert>
#include <math.h>
#include <pthread.h>
#include "printLock.h"
#include "parseInput.h"
#include "LinearSmithWaterman.h"
#include "LinearNeedlemanWunsch.h"
#include "AffineNeedlemanWunsch.h"
#include "timing.h"

// DEBUG printing for multithreading
//#define DEBUG

// MULTHREADED DEFINE (USE_THREADS) IS IN PRINTLOCK.H
// REMEMBER THAT MULTITHREADING OUTPUT ORDERING IS NONDETERMINISTIC -> USE reorderOutput.py [input file] [output file]
// Multithreaded parameters 
#define PAIRS_PER_THREAD 20
#define THREADS_PER_BATCH 20

// Choose what alignment we are doing (ONLY CHOOSE ONE)
#define LSW_ENABLE
//#define LNW_ENABLE
//#define ANW_ENABLE


using std::strcmp;

struct thread_arg{
    int threadID;
    int batchID;
    int pairNum;
    int numPairs;
    char* sequences;
    seqPair* sequenceIdxs;
    int matchWeight;
    int mismatchWeight;
    int gapOpenWeight;
    int gapExtendWeight;
};

void* threadComputeLSW(void *tmp){
	struct thread_arg* threadData = (struct thread_arg*)tmp;
    int threadID            = threadData->threadID;
    int batchID             = threadData->batchID;
    int pairNum             = threadData->pairNum;
    int numPairs            = threadData->numPairs;
    char* sequences         = threadData->sequences;
	seqPair* sequenceIdxs   = threadData->sequenceIdxs;
    int matchWeight         = threadData->matchWeight;
    int mismatchWeight      = threadData->mismatchWeight;
    int gapWeight           = threadData->gapOpenWeight;

    #ifdef DEBUG
        printLock();
        printf("Thread ID: %d Batch ID: %d\n", threadID, batchID);
        fflush(stdout);
        printUnlock();
    #endif

    for(int i = pairNum; i < pairNum + numPairs; i++){
        LinearSmithWaterman LSW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], i, matchWeight, mismatchWeight, gapWeight);
        LSW.align();
    }
}

void* threadComputeLNW(void *tmp){
	struct thread_arg* threadData = (struct thread_arg*)tmp;
    int threadID            = threadData->threadID;
    int batchID             = threadData->batchID;
    int pairNum             = threadData->pairNum;
    int numPairs            = threadData->numPairs;
    char* sequences         = threadData->sequences;
	seqPair* sequenceIdxs   = threadData->sequenceIdxs;
    int matchWeight         = threadData->matchWeight;
    int mismatchWeight      = threadData->mismatchWeight;
    int gapWeight           = threadData->gapOpenWeight;

    #ifdef DEBUG
        printLock();
        printf("Thread ID: %d Batch ID: %d\n", threadID, batchID);
        fflush(stdout);
        printUnlock();
    #endif

    for(int i = pairNum; i < pairNum + numPairs; i++){
        LinearNeedlemanWunsch LNW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], i, matchWeight, mismatchWeight, gapWeight);
        LNW.align();
    }
}

void* threadComputeANW(void *tmp){
	struct thread_arg* threadData = (struct thread_arg*)tmp;
    int threadID            = threadData->threadID;
    int batchID             = threadData->batchID;
    int pairNum             = threadData->pairNum;
    int numPairs            = threadData->numPairs;
    char* sequences         = threadData->sequences;
	seqPair* sequenceIdxs   = threadData->sequenceIdxs;
    int matchWeight         = threadData->matchWeight;
    int mismatchWeight      = threadData->mismatchWeight;
    int gapOpenWeight       = threadData->gapOpenWeight;
    int gapExtendWeight     = threadData->gapExtendWeight;

    #ifdef DEBUG
        printLock();
        printf("Thread ID: %d Batch ID: %d\n", threadID, batchID);
        fflush(stdout);
        printUnlock();
    #endif

    for(int i = pairNum; i < pairNum + numPairs; i++){
        AffineNeedlemanWunsch ANW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], i, matchWeight, mismatchWeight, gapOpenWeight, gapExtendWeight);
        ANW.align();
    }
}

int main(int argc, char *argv[]){

    // Check that YOU use it correctly
    if (argc < 2) {
		fprintf(stderr, "usage: main -pairs <InSeqFile> -match <matchWeight> -mismatch <mismatchWeight> -gap <gapWeight> \n");
		exit(EXIT_FAILURE);
    }
	
    // Get filename from args
    char *pairFileName;
    int matchWeight     = 3;
    int mismatchWeight  = -1;
    // int gapWeight   = -2;
    int gapOpenWeight   = -4;
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
    // if(argc > 7 && strcmp(argv[7], "-gap") == 0) {
    //     gapWeight = atoi(argv[8]);
    // }
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

    // Print out the first byte of each section to make sure I'm not crazy
    // std::cout << "Printing sequence info" << std::endl;
    // printParsedFile(fileInfo.numPairs, sequenceIdxs, sequences);

    // Setup threads
    uint64_t start_time = start_timer();
    printf("Pair # | Score\n");
    #ifdef USE_THREADS
        int numPairsPerThread   = PAIRS_PER_THREAD;
        int threadsPerBatch     = THREADS_PER_BATCH;
        int numBatches          = ceil((fileInfo.numPairs/numPairsPerThread)/threadsPerBatch);
        int batchID;
        int threadID;
        for(batchID = 0; batchID < numBatches; batchID++){
            int numThreadsComputing = 0;
            pthread_t           *threads    = (pthread_t*)malloc(sizeof(pthread_t)*threadsPerBatch);
            struct thread_arg   *threadArgs = (thread_arg*)malloc(sizeof(thread_arg)*threadsPerBatch);
            for(threadID = 0; threadID < threadsPerBatch; threadID++){
                struct thread_arg threadArg;
                threadArg.threadID  = threadID;
                threadArg.batchID   = batchID;
                threadArg.pairNum   = batchID*threadsPerBatch*numPairsPerThread + threadID*numPairsPerThread;
                // Dont create a thread if we exceed the number of pairs
                if(threadArg.pairNum > fileInfo.numPairs - 1){
                    break;
                }
                threadArg.numPairs        = numPairsPerThread;
                threadArg.sequences       = sequences;
                threadArg.sequenceIdxs    = sequenceIdxs;
                threadArg.matchWeight     = matchWeight;
                threadArg.mismatchWeight  = mismatchWeight;
                threadArg.gapOpenWeight   = gapOpenWeight;
                threadArg.gapExtendWeight = gapExtendWeight;
                threadArgs[threadID]      = threadArg;
                
                #ifdef DEBUG
                    printLock();
                    printf("Attempting to create thread:%d\n", threadID);
                    fflush(stdout);
                    printUnlock();
                #endif
                
                int ret;
                #ifdef LSW_ENABLE
                    ret = pthread_create(&threads[threadID], NULL, threadComputeLSW, (void*) &threadArgs[threadID]);
                #endif
                #ifdef LNW_ENABLE
                    ret = pthread_create(&threads[threadID], NULL, threadComputeLNW, (void*) &threadArgs[threadID]);
                #endif
                #ifdef ANW_ENABLE
                    ret = pthread_create(&threads[threadID], NULL, threadComputeANW, (void*) &threadArgs[threadID]);
                #endif
                assert(ret == 0);
                numThreadsComputing++;
            }

            // Wait for each thread to terminate
            for(threadID = 0; threadID < numThreadsComputing; threadID++){
                #ifdef DEBUG
                    printLock();
                    printf("Attempting to wait for thread thread:%d\n", threadID);
                    fflush(stdout);
                    printUnlock();
                #endif
                int ret;
                ret = pthread_join(threads[threadID], NULL);
                assert(ret == 0);
            }

            // Free thread arrays
            free(threads);
            free(threadArgs);

        }

    #else
        // Score some matrices babyyyyyyyyyyy
        #ifdef LSW_ENABLE
            for(int i = 0; i < fileInfo.numPairs; i++){
                LinearSmithWaterman LSW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], i, matchWeight, mismatchWeight, gapOpenWeight);
                LSW.align();
            }
        #endif
        #ifdef LNW_ENABLE
            for(int i = 0; i < fileInfo.numPairs; i++){
                LinearNeedlemanWunsch LNW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], i, matchWeight, mismatchWeight, gapOpenWeight);
                LNW.align();
            }
        #endif
        #ifdef ANW_ENABLE
            for(int i = 0; i < fileInfo.numPairs/; i++){
                AffineNeedlemanWunsch ANW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], i, matchWeight, mismatchWeight, gapOpenWeight, gapExtendWeight);
                ANW.align();
            }
        #endif
    #endif

    uint64_t elapsed_time = get_elapsed_time();
    printf("Elapsed time (usec): %lld\n", elapsed_time);

    // Free data arrays
    printf("Cleaning up\n");
    cleanupParsedFile(sequenceIdxs, sequences);
}