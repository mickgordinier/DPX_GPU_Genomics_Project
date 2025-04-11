#include <cstring>
#include "parseInput.h"
#include "LinearSmithWaterman.h"
#include "LinearNeedlemanWunsch.h"

using std::strcmp;

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

    // Print out the first byte of each section to make sure I'm not crazy
    // std::cout << "Printing sequence info" << std::endl;
    // printParsedFile(numPairs, sequenceIdxs, sequences);

    // Score some matrices babyyyyyyyyyyy
    printf("Pair # | Score\n");
    // LinearSmithWaterman LSW("GAGTACTCAACACCAACATTGATGGGCAATGGAAAATAGCCTTCGCCATCACACCATTAAGGGTGA", "GAATACTCAACAGCAACATCAACGGGCAGCAGAAAATAGGCTTTGCCATCACTGCCATTAAGGATGTGGG", 3, -1, -2);
    // for(int i = 0; i < numPairs/1000; i++){
    //     printf("%d | ", i);
    //     LinearSmithWaterman LSW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], matchWeight, mismatchWeight, gapWeight);
    //     LSW.align();
    // }

    // LinearNeedlemanWunsch LNW("GAGTACTCAACACCAACATTGATGGGCAATGGAAAATAGCCTTCGCCATCACACCATTAAGGGTGA", "GAATACTCAACAGCAACATCAACGGGCAGCAGAAAATAGGCTTTGCCATCACTGCCATTAAGGATGTGGG", 3, -1, -2);
    for(int i = 0; i < numPairs/1000; i++){
        printf("%d |", i);
        LinearNeedlemanWunsch LNW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], matchWeight, mismatchWeight, gapWeight);
        LNW.align();
    }

    // Free data arrays
    printf("Cleaning up\n");
    cleanupParsedFile(sequenceIdxs, sequences);
}