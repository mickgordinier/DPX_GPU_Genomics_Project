#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <algorithm>
#include <cstdint>

#define PRINT_PARSED_PAIRS

struct inputInfo{
    size_t numPairs;
    size_t numBytes;
    size_t numCells;
    size_t minReferenceLength;
    size_t minQueryLength;
    size_t maxReferenceLength;
    size_t maxQueryLength;
    double avgReferenceLength;
    double avgQueryLength;

};

struct seqPair {
    // size_t scoreSeedIdx;
    // size_t scoreSeedSize;
    int referenceIdx;
    int referenceSize;
    int queryIdx;
    int querySize;
};

inputInfo parseInput(const char* pairFileName, seqPair* &sequence_indices, char* &sequences);

void printParsedFile(const size_t numPairs, const seqPair* sequence_indices, const char* sequences);

void cleanupParsedFile(seqPair* sequence_indices, char* sequences);