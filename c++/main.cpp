#include <iostream>
#include <cstring>
#include "parseInput.h"
#include "LinearSmithWaterman.h"
#include "LinearNeedlemanWunsch.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::getline;
using std::strcmp;



int main(int argc, char *argv[]){

    // Check that YOU use it correctly
    if (argc < 2) {
		fprintf(stderr, "usage: main -pairs <InSeqFile>\n");
		exit(EXIT_FAILURE);
	}
	
    // Get filename from args
    string pairFileName;
	if(strcmp(argv[1], "-pairs") == 0) {
		pairFileName = argv[2];
    }

    // Parse input file
    std::cout << "Parsing input file: " << pairFileName << std::endl;
    size_t numPairs;
    seqPair* sequenceIdxs;
    char* sequences;
	numPairs = parseInput(pairFileName.c_str(), sequenceIdxs, sequences);

    // Print out the first byte of each section to make sure I'm not crazy
    // std::cout << "Printing sequence info" << std::endl;
    // printParsedFile(numPairs, sequenceIdxs, sequences);

    // Score some matrices babyyyyyyyyyyy
    std::cout << "Pair # | Score\n";
    // LinearSmithWaterman LSW("GAGTACTCAACACCAACATTGATGGGCAATGGAAAATAGCCTTCGCCATCACACCATTAAGGGTGA", "GAATACTCAACAGCAACATCAACGGGCAGCAGAAAATAGGCTTTGCCATCACTGCCATTAAGGATGTGGG", 3, -1, -2);
    for(int i = 0; i < numPairs/1000; i++){
        std::cout << i << " | ";
        LinearSmithWaterman LSW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], 3, -1, -2);
        LSW.align();
    }

    // LinearNeedlemanWunsch LNW("GAGTACTCAACACCAACATTGATGGGCAATGGAAAATAGCCTTCGCCATCACACCATTAAGGGTGA", "GAATACTCAACAGCAACATCAACGGGCAGCAGAAAATAGGCTTTGCCATCACTGCCATTAAGGATGTGGG", 3, -1, -2);
    for(int i = 0; i < numPairs/1000; i++){
        std::cout << i << " | ";
        LinearNeedlemanWunsch LNW(&sequences[sequenceIdxs[i].referenceIdx], &sequences[sequenceIdxs[i].queryIdx], 3, -1, -2);
        LNW.align();
    }

    // Free data arrays
    std::cout << "Cleanup data arrays" << std::endl;
    cleanupParsedFile(sequenceIdxs, sequences);
}