/*
Description: Class for performing linear-gap penalty global Needleman-Wunsch alsignment
*/
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include "SequenceAligner.h"
#include "debug.h"
#include "printLock.h"

class AffineNeedlemanWunsch: public SequenceAligner {
    
    private:
    
        // Weights for sequence scoring
        int matchWeight;
        int mismatchWeight;
        int gapOpenWeight;
        int gapExtendWeight;

        // Scoring matrix
        std::vector<std::vector<int> > scoringMemo;
        std::vector<std::vector<int> > queryInsertionMemo;
        std::vector<std::vector<int> > queryDeletionMemo;
    
        // Enumerated predecessor directions for backtracking
        enum directionMain {
            NONE_MAIN,             // For initialization purposes
            MATCH,
            MISMATCH,
            QUERY_INSERTION,
            QUERY_DELETION
        }; 

        // If in one of the other matrices besides the main scoring matrix
        // Element can be computed from either:
        //   i. Opening a gap from the previous element
        //  ii. Extending the gap from
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

        // Matrix of predecessor directions for every corresponding cell in memo
        std::vector<std::vector<directionMain> > scoringBacktrack;
        std::vector<std::vector<directionIndel> > queryInsertionBacktrack;
        std::vector<std::vector<directionIndel> > queryDeletionBacktrack;

    public:
    
        AffineNeedlemanWunsch(const std::string inputReference, const std::string inputQuery, const int pairNum,
                                const int matchWeight, const int mismatchWeight, 
                                const int gapOpenWeight, const int gapExtendWeight) : 
            SequenceAligner(inputReference, inputQuery, pairNum), 
            matchWeight(matchWeight), 
            mismatchWeight(mismatchWeight), 
            gapOpenWeight(gapOpenWeight),
            gapExtendWeight(gapExtendWeight) {}

        void init_matrix(); // Initialies matrices

        void print_matrix(); // Prints current scoring matrix

        void score_matrix(); // Scores the matrix

        void backtrack(); // Backtracks to find highest scoring sequences

        void align(); // Performs entire alignment, from initialization to result printing

        void print_results(); // Prints results
};