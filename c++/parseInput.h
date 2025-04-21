#include <cstdio>
#include <cstdlib>
#include <cstddef>

#define PRINT_PARSED_PAIRS

struct inputInfo{
    size_t numPairs;
    size_t numBytes;
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