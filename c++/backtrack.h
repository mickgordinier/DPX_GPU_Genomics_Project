#include <iostream>
#include <string>

// Enumerated predecessor directions for backtracking
enum direction {
    NONE,             // For initialization purposes
    MATCH,
    MISMATCH,
    QUERY_INSERTION,
    QUERY_DELETION
};  

void printMatrix(const int *memo, const int referenceLength, const int queryLength);

void backtrackNW(const direction *backtrackMemo, const char *referenceString, const int referenceLength, const char *queryString,  const int queryLength);