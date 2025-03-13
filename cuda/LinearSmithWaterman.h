#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include "SequenceAligner.h"
#include "debug.h"

class LinearSmithWaterman: public SequenceAligner{
    private:
    // [Used for scoring]
    int match_weight;
    int mismatch_weight;
    int gap_weight;
    std::vector<std::vector<int> > memo; // Matrix containing all scores
    
    // [Used for backtracking]
    int max_score;

    // Enumerated predecessor directions for backtracking
    enum direction {
        NONE,
        CORNER_MATCH,
        CORNER_MISMATCH,
        LEFT_GAP,
        UPPER_GAP
    }; 

    // Matrix of predecessor directions for every corresponding cell in memo
    std::vector<std::vector<direction> > backtrack_memo;

    // Struct used for keeping track of backtracking info
    struct backtrack_info {
        std::string query_sequence;
        std::string pair_relation;
        std::string reference_sequence;
        size_t row_idx;
        size_t col_idx;
    };

    // Queue that we add to and pull from to backtrack cells
    std::deque<backtrack_info> backtrack_queue;

    // Vector of all of our final subsequences
    std::vector<backtrack_info> results;


    public:
    LinearSmithWaterman(const std::string input_reference, const std::string input_query, const int match_weight, const int mismatch_weight, const int gap_weight) : 
        SequenceAligner(input_reference, input_query), 
        match_weight(match_weight), 
        mismatch_weight(mismatch_weight), 
        gap_weight(gap_weight) {
            max_score = 0;
        };

    void init_matrix(); // Initialies matrices

    void print_matrix(); // Prints current scoring matrix

    void score_matrix(); // Scores the matrix

    void backtrack(); // Backtracks to find highest scoring sequences

    void align(); // Performs entire alignment, from initialization to result printing

    void print_results(); // Prints results

};