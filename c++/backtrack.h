#include <iostream>
#include <string>
#include  "printLock.h"

// Enumerated predecessor directions for backtracking
// enum direction {
//     NONE,             // For initialization purposes
//     MATCH,
//     MISMATCH,
//     QUERY_INSERTION,
//     QUERY_DELETION
// };  

enum directionMain {
    NONE_MAIN,             // For initialization purposes
    MATCH,
    MISMATCH,
    QUERY_INSERTION,
    QUERY_DELETION
};

// For Affine NW Backtracking
enum directionIndel {
    NONE_INDEL,
    GAP_OPEN,
    GAP_EXTEND
};

enum currentMatrixPosition {
    SCORING,
    INSERTION,
    DELETION
};

void printMatrix(const int *memo, const int referenceLength, const int queryLength);

void printBacktrackMatrix(const directionMain *memo, const int referenceLength, const int queryLength);

void backtrackNW(const directionMain *backtrackMemo, const char *referenceString, const int referenceLength, const char *queryString,  const int queryLength);

void backtrackMultiNW(const directionMain *backtrackMemo, const char *referenceString, const int referenceLength, const char *queryString, const int queryLength, const int pairNum, const int score);

void backtrackANW(
    const directionMain *scoringBacktrack, 
    const directionIndel *queryInsertionBacktrack,
    const directionIndel *queryDeletionBacktrack,
    const char *referenceString, const int referenceLength, 
    const char *queryString,  const int queryLength);