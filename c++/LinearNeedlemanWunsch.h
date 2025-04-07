/*
Description: Class for performing linear-gap penalty global Needleman-Wunsch alsignment
*/
#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include "SequenceAligner.h"
#include "debug.h"

class LinearNeedlemanWunsch: public SequenceAligner {
    
    private:
    
        // Weights for sequence scoring
        int match_weight;
        int mismatch_weight;
        int gap_weight;

        // Scoring matrix
        std::vector<std::vector<int> > memo;
    
        // Enumerated predecessor directions for backtracking
        enum direction {
            NONE,             // For initialization purposes
            MATCH,
            MISMATCH,
            QUERY_INSERTION,
            QUERY_DELETION
        }; 

        // Matrix of predecessor directions for every corresponding cell in memo
        std::vector<std::vector<direction> > backtrack_memo;

        // Vector of all of our final subsequences
        // std::vector<backtrack_info> results;


    public:
    
        LinearNeedlemanWunsch(const std::string input_reference, const std::string input_query, 
                                const int match_weight, const int mismatch_weight, const int gap_weight) : 
            SequenceAligner(input_reference, input_query), 
            match_weight(match_weight), 
            mismatch_weight(mismatch_weight), 
            gap_weight(gap_weight) {}

        void init_matrix(); // Initialies matrices

        void print_matrix(); // Prints current scoring matrix

        void score_matrix(); // Scores the matrix

        void backtrack(); // Backtracks to find highest scoring sequences

        void align(); // Performs entire alignment, from initialization to result printing

        void print_results(); // Prints results
};