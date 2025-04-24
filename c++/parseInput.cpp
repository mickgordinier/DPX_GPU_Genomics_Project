#include "parseInput.h"

#define PARSE_SCORE_SEED 0
#define PARSE_REFERENCE 1
#define PARSE_QUERY 2

#define INPUT_CAP 10000000

inputInfo parseInput(const char* pairFileName, seqPair* &sequenceIdxs, char* &sequences){
    // Open file
    FILE *pairFile = fopen(pairFileName, "r");	
    if (pairFile == NULL) {
        fprintf(stderr, "Could not open file: %s\n", pairFileName);
        exit(1);
    }

    // Get number of lines in file and setup indices into data array
    const int bufSize = 1024 * 1024;
    char* buffer = (char*)malloc(bufSize * sizeof(char));
    size_t numLines = 0;
    size_t numBytes = 0;
    size_t n;
    while (n = fread(buffer, sizeof(char), bufSize, pairFile)) {
        if(n >= 0){
            numBytes += n;
        } else {
            fprintf(stderr, "Error occurred while reading file: %d\n", n);
            exit(1);
        }
        for (int i = 0; i < n; i++) {
            if (buffer[i] == '\n') {
                numLines++;
            }
        }
    }

    // Calculate number of pairs
    if(numLines % 3 != 0){
        fprintf(stderr, "Number of lines not a multiple of 3: %s\n", numLines);
        exit(1);
    }
    size_t numPairs = numLines/3;

    // Free buffer
    free(buffer);

    // Allocate memory
    sequenceIdxs = (seqPair*)malloc(numPairs * sizeof(seqPair));
    sequences        = (char*)malloc(numBytes * sizeof(char));

    // Reset file pointer and read in the entire file data
    fseek(pairFile, 0L, SEEK_SET);
    size_t bytesRead = 0;
    while(n = fread(sequences + bytesRead, sizeof(char), numBytes - bytesRead, pairFile)){
        bytesRead += n;
        if(!(n >= 0)){
            fprintf(stderr, "Error occurred while reading file second time: %d\n", n);
            exit(1);
        }
    }
    if(bytesRead != numBytes){
        fprintf(stderr, "Did not read all bytes: %d read, %d total\n", bytesRead, numBytes);
        exit(1);
    }
    
    inputInfo retVal;
    retVal.numCells             = 0;
    retVal.minReferenceLength   = SIZE_MAX;
    retVal.minQueryLength       = SIZE_MAX;
    retVal.maxReferenceLength   = 0;
    retVal.maxQueryLength       = 0;
    retVal.avgReferenceLength   = 0;
    retVal.avgQueryLength       = 0;


    // Read all of the bytes into our data struct
    int parseMode = PARSE_SCORE_SEED;
    size_t sequenceIdx = 0;
    // sequenceIdxs[sequenceIdx].scoreSeedIdx = 0;
    for(size_t i = 0; i < numBytes; i++){
        if (sequences[i] == '\n') {
            sequences[i] = '\0'; // replace with null terminator
            if(parseMode == PARSE_SCORE_SEED){
                // sequenceIdxs[sequenceIdx].scoreSeedSize = i - sequenceIdxs[sequenceIdx].scoreSeedIdx;
                sequenceIdxs[sequenceIdx].referenceIdx = i + 1;
                parseMode = PARSE_REFERENCE;
            } else if (parseMode == PARSE_REFERENCE){
                sequenceIdxs[sequenceIdx].referenceSize = i - sequenceIdxs[sequenceIdx].referenceIdx;
                retVal.avgReferenceLength += sequenceIdxs[sequenceIdx].referenceSize;
                retVal.maxReferenceLength = std::max(retVal.maxReferenceLength, (size_t)sequenceIdxs[sequenceIdx].referenceSize);
                retVal.minReferenceLength = std::min(retVal.minReferenceLength, (size_t)sequenceIdxs[sequenceIdx].referenceSize);
                sequenceIdxs[sequenceIdx].queryIdx = i + 1;
                parseMode = PARSE_QUERY;
            } else if (parseMode == PARSE_QUERY){
                sequenceIdxs[sequenceIdx].querySize = i - sequenceIdxs[sequenceIdx].queryIdx;
                retVal.avgQueryLength += sequenceIdxs[sequenceIdx].querySize;
                retVal.maxQueryLength = std::max(retVal.maxQueryLength, (size_t)sequenceIdxs[sequenceIdx].querySize);
                retVal.minQueryLength = std::min(retVal.minQueryLength, (size_t)sequenceIdxs[sequenceIdx].querySize);
                retVal.numCells += (sequenceIdxs[sequenceIdx].referenceSize * sequenceIdxs[sequenceIdx].querySize);
                sequenceIdx += 1;
                if(sequenceIdx == INPUT_CAP){
                    retVal.numPairs = INPUT_CAP;
                    break;
                }
                parseMode = PARSE_SCORE_SEED;
                // sequenceIdxs[sequenceIdx].scoreSeedIdx = i + 1;
            } else {
                fprintf(stderr, "Parse mode tweaking: %d\n", n);
                exit(1);
            }
        }
    }
    retVal.numBytes = numBytes;
    retVal.avgReferenceLength = retVal.avgReferenceLength / retVal.numPairs;
    retVal.avgQueryLength = retVal.avgQueryLength / retVal.numPairs;
    
    return retVal;
}

void printParsedFile(const size_t numPairs, const seqPair* sequenceIdxs, const char* sequences){
// Print out the first byte of each section to make sure I'm not crazy
    size_t sequenceIdx = 0;
    while(true){
        printf("Pair: %d Reference: %c, Reference Size: %d, Query: %c, Query Size: %d\n", 
        sequenceIdx,
        // sequences[sequenceIdxs[sequenceIdx].scoreSeedIdx],
        // sequenceIdxs[sequenceIdx].scoreSeedSize,
        sequences[sequenceIdxs[sequenceIdx].referenceIdx],
        sequenceIdxs[sequenceIdx].referenceSize,
        sequences[sequenceIdxs[sequenceIdx].queryIdx],
        sequenceIdxs[sequenceIdx].querySize);
        sequenceIdx += 1;
        if(sequenceIdx == numPairs){
            break;
        }
    }
}

void cleanupParsedFile(seqPair* sequenceIdxs, char* sequences){
    free(sequenceIdxs);
    free(sequences);
}