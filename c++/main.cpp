#include <cstring>
#include <cassert>
#include <pthread.h>
#include "printLock.h"
#include "parseInput.h"
#include "LinearSmithWaterman.h"
#include "LinearNeedlemanWunsch.h"
#include "AffineNeedlemanWunsch.h"
#include "timing.h"

//#define DEBUG

using std::strcmp;

struct thread_arg{
    int threadID;
    int batchID;
    int pairNum;
    char* reference;
    char* query;
    int matchWeight;
    int mismatchWeight;
    int gapOpenWeight;
    int gapExtendweight;
};

void* threadComputeLNW(void *tmp){
	struct thread_arg* threadData = (struct thread_arg*)tmp;
    int threadID            = threadData->threadID;
    int batchID             = threadData->batchID;
    int pairNum             = threadData->pairNum;
    char* reference         = threadData->reference;
	char* query             = threadData->query;
    int matchWeight         = threadData->matchWeight;
    int mismatchWeight      = threadData->mismatchWeight;
    int gapWeight           = threadData->gapOpenWeight;

    #ifdef DEBUG
        printLock();
        printf("Thread ID: %d Batch ID: %d\n", threadID, batchID);
        fflush(stdout);
        printUnlock();
    #endif

    LinearNeedlemanWunsch LNW(reference, query, pairNum, matchWeight, mismatchWeight, gapWeight);
    LNW.align();
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
        int threadsPerBatch = 100;
        int numBatches      = fileInfo.numPairs/100;
        int batchID;
        int threadID;
        for(batchID = 0; batchID < numBatches; batchID++){
            pthread_t           *threads    = (pthread_t*)malloc(sizeof(pthread_t)*threadsPerBatch);
            struct thread_arg   *threadArgs = (thread_arg*)malloc(sizeof(thread_arg)*threadsPerBatch);
            for(threadID = 0; threadID < threadsPerBatch; threadID++){
                struct thread_arg threadArg;
                threadArg.threadID  = threadID;
                threadArg.batchID   = batchID;
                threadArg.pairNum   = batchID*threadsPerBatch + threadID;
                threadArg.reference = &sequences[sequenceIdxs[batchID*threadsPerBatch + threadID].referenceIdx];
                threadArg.query     = &sequences[sequenceIdxs[batchID*threadsPerBatch + threadID].queryIdx];
                threadArg.matchWeight       = matchWeight;
                threadArg.mismatchWeight    = mismatchWeight;
                threadArg.gapOpenWeight     = gapOpenWeight;
                threadArgs[threadID]        = threadArg;
                
                #ifdef DEBUG
                    printLock();
                    printf("Attempting to create thread:%d\n", threadID);
                    fflush(stdout);
                    printUnlock();
                #endif
                
                int ret;
                ret = pthread_create(&threads[threadID], NULL, threadComputeLNW, (void*) &threadArgs[threadID]);
                assert(ret == 0);
            }

            // Wait for each thread to terminate
            for(threadID = 0; threadID < threadsPerBatch; threadID++){
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
        // LinearSmithWaterman LSW("GAGTACTCAACACCAACATTGATGGGCAATGGAAAATAGCCTTCGCCATCACACCATTAAGGGTGA", "GAATACTCAACAGCAACATCAACGGGCAGCAGAAAATAGGCTTTGCCATCACTGCCATTAAGGATGTGGG", 3, -1, -2);
        // for(int i = 0; i < fileInfo.numPairs; i++){
        //     printf("%d | ", i);
        //     LinearSmithWaterman LSW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], matchWeight, mismatchWeight, gapOpenWeight);
        //     LSW.align();
        // }

        for(int i = 0; i < fileInfo.numPairs; i++){
            //printf("%d | ", i);
            LinearNeedlemanWunsch LNW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], i, matchWeight, mismatchWeight, gapOpenWeight);
            LNW.align();
        }

        // for(int i = 0; i < fileInfo.numPairs/100; i++){
        //     printf("%d | ", i);
        //     AffineNeedlemanWunsch ANW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], matchWeight, mismatchWeight, gapOpenWeight, gapExtendWeight);
        //     ANW.align();
        // }
    #endif

    uint64_t elapsed_time = get_elapsed_time();
    printf("Elapsed time (usec): %lld\n", elapsed_time);

    // Free data arrays
    printf("Cleaning up\n");
    cleanupParsedFile(sequenceIdxs, sequences);
}