#include <iostream>
#include <iomanip>
#include <vector>
#include <deque>
#include "SequenceAligner.h"
#include "debug.h"

class LinearSmithWaterman: public SequenceAligner{
    private:
    // Used for scoring
    int match_weight;
    int mismatch_weight;
    int gap_weight;
    std::vector<std::vector<int> > memo;

    // Used for backtracking
    enum direction {
        NONE,
        CORNER_MATCH,
        CORNER_MISMATCH,
        LEFT_GAP,
        UPPER_GAP
    };

    std::vector<std::vector<direction> > backtrack_memo;

    struct backtrack_info {
        std::string query_sequence;
        std::string pair_relation;
        std::string reference_sequence;
        size_t row_idx;
        size_t col_idx;
        int score;
    };

    std::deque<backtrack_info> backtrack_queue;

    std::vector<backtrack_info> results;


    public:
    LinearSmithWaterman(const std::string input_reference, const std::string input_query, const int match_weight, const int mismatch_weight, const int gap_weight) : 
        SequenceAligner(input_reference, input_query), 
        match_weight(match_weight), 
        mismatch_weight(mismatch_weight), 
        gap_weight(gap_weight) {};

    void init_matrix();

    void print_matrix();

    void score_matrix();

    void backtrack();

    void align();


    void print_results();

};