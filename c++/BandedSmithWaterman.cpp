#include "BandedSmithWaterman.h"

using std::cout;
using std::vector;
using std::max;
using std::string;

void BandedSmithWaterman::init_matrix(){
    #ifdef PRINT_EXTRA
        cout << "[Initializing Matrix...]\n";
    #endif
    size_t num_rows = query_str.size();
    size_t num_cols = reference_str.size();

    // Initialize memo and backtrack with only band width
    vector<int> temp_memo(num_cols + 1, 0);
    memo.resize(num_rows + 1, temp_memo);

    vector<direction> temp_backtrack(num_cols, NONE);
    backtrack_memo.resize(num_rows, temp_backtrack);
}

void BandedSmithWaterman::print_matrix(){
    size_t num_rows = query_str.size();
    size_t num_cols = reference_str.size();

    // Print parameters
    cout << "Reference: " << reference_str << " Size: " << num_cols << "\n";
    cout << "Query: " << query_str << " Size: " << num_rows << "\n";

    cout << "Matrix Dim: [ " << num_rows + 1 << " x " << num_cols + 1 << " ]\n";

    // Set the spacing for each element of the matrix
    int min_width = 2;

    // Print out reference str
    cout << "    " << std::setw(min_width) << " " << " ";
    for(size_t i = 0; i < reference_str.size(); i++){
        cout << "  " << std::setw(min_width) << reference_str[i] << " ";
    }
    cout << "\n";

    // Print rest of matrix
    for(size_t i = 0; i < memo.size(); i++){
        if(i == 0){
            cout << "  ";
        }
        else {
            cout << query_str[i-1] << " ";
        }
        cout << "[";
        for(size_t j = 0; j < memo[i].size(); j++){
            if(j != memo[i].size() - 1){
                cout << " " << std::setw(min_width) << memo[i][j] << ", ";
            }
            else {
                cout << " " << std::setw(min_width) << memo[i][j];
            }
        } //end for
        cout << "],\n";
    } //end for

    cout << std::endl;

}

void BandedSmithWaterman::score_matrix() {
    #ifdef PRINT_EXTRA
        cout << "[Scoring (banded)...]\n";
    #endif

    for(size_t row_idx = 1; row_idx < memo.size(); row_idx++) {
        // Only compute where abs(row - col) <= 1
        
        size_t start_col = 0;
        size_t end_col = 0;

        //we are doing computation from one column left of the diagonal
        if (row_idx >= 1) {
            start_col = row_idx - band_width;
        } else {
            start_col = band_width;
        }
        
        end_col = std::min(row_idx + 1, memo[0].size() - 1); //we are trying to do computation one column to the right of the diagonal, memo[0].size() - 1 indicates this last column

        for(size_t col_idx = start_col; col_idx <= end_col; col_idx++) {
            int upper_score;
            int left_score;
            int corner_score;

            for (size_t col_idx = start_col; col_idx <= end_col; col_idx++) {
                int upper_score = memo[row_idx - 1][col_idx] + gap_weight;
                int left_score = memo[row_idx][col_idx - 1] + gap_weight;
                int corner_score;
                direction corner_direction;

                if (query_str[row_idx - 1] == reference_str[col_idx - 1]) {
                    corner_score = memo[row_idx - 1][col_idx - 1] + match_weight;
                    corner_direction = CORNER_MATCH;
                } else {
                    corner_score = memo[row_idx - 1][col_idx - 1] + mismatch_weight;
                    corner_direction = CORNER_MISMATCH;
                }

                int temp_greatest = max(upper_score, max(left_score, corner_score));

                memo[row_idx][col_idx] = max(0, temp_greatest);

                if (temp_greatest == corner_score) {
                    backtrack_memo[row_idx - 1][col_idx - 1] = corner_direction;
                } else if (temp_greatest == left_score) {
                    backtrack_memo[row_idx - 1][col_idx - 1] = LEFT_GAP;
                } else {
                    backtrack_memo[row_idx - 1][col_idx - 1] = UPPER_GAP;
                }
            }
        }
    }
}

void BandedSmithWaterman::backtrack(){
    for(size_t row_idx = memo.size() - 1; row_idx > 0; row_idx--) {
        size_t min_col = (row_idx > band_width) ? row_idx - band_width : 1;
        size_t max_col = std::min(row_idx + band_width, memo[0].size() - 1);
        
        for(size_t col_idx = memo[0].size() - 1; col_idx > 0; col_idx--) {
            int new_max_score = max(max_score, memo[row_idx][col_idx]);
            if (new_max_score > max_score) {
                backtrack_queue.clear();
                max_score = new_max_score;
            }

            if (memo[row_idx][col_idx] == max_score) {
                backtrack_info info = {"", "", "", row_idx, col_idx};
                backtrack_queue.push_back(info);
            }

            if (col_idx == min_col){
                break;
            }
        }
    }
    // Now do the actual backtracking
    while(!backtrack_queue.empty()){
        // Remove a cell from our queue
        backtrack_info current_cell = backtrack_queue.front();
        backtrack_queue.pop_front();

        backtrack_info next_cell;
        // Determine the current cell's predecessor
        switch (backtrack_memo[current_cell.row_idx-1][current_cell.col_idx-1]) {
            case CORNER_MATCH:   
                next_cell.col_idx = current_cell.col_idx-1;
                next_cell.row_idx = current_cell.row_idx-1;
                next_cell.reference_sequence = reference_str[current_cell.col_idx-1];
                next_cell.pair_relation = "*";
                next_cell.query_sequence = query_str[current_cell.row_idx-1];
                break;
            // end if match
            
            case CORNER_MISMATCH:
                next_cell.col_idx = current_cell.col_idx-1;
                next_cell.row_idx = current_cell.row_idx-1;
                next_cell.reference_sequence = reference_str[current_cell.col_idx-1];
                next_cell.pair_relation = "|";
                next_cell.query_sequence = query_str[current_cell.row_idx-1];
                break;
            // end if mismatch
            
            case LEFT_GAP:
                next_cell.col_idx = current_cell.col_idx-1;
                next_cell.row_idx = current_cell.row_idx;
                next_cell.reference_sequence = reference_str[current_cell.col_idx-1];
                next_cell.pair_relation = " ";
                next_cell.query_sequence = "_";
                break;
            // end if left gap
            
            case UPPER_GAP:
                next_cell.col_idx = current_cell.col_idx;
                next_cell.row_idx = current_cell.row_idx-1;
                next_cell.reference_sequence = "_";
                next_cell.pair_relation = " ";
                next_cell.query_sequence = query_str[current_cell.row_idx-1];
                break;
            default:
                break;
            // end if upper gap
        } // end switch

        // Update our accumulated sequences with the current cell's reference + query index
        //DEBUG_PRINT("backtrack", "Reference: " << next_cell.reference_sequence);
        //DEBUG_PRINT("backtrack", "Relation: " << next_cell.pair_relation);
        //DEBUG_PRINT("backtrack", "Query: " << next_cell.query_sequence);
        next_cell.reference_sequence = next_cell.reference_sequence + current_cell.reference_sequence;
        next_cell.pair_relation = next_cell.pair_relation + current_cell.pair_relation;
        next_cell.query_sequence = next_cell.query_sequence + current_cell.query_sequence;
        //DEBUG_PRINT("backtrack", "Cumulative Reference: " << next_cell.reference_sequence);
        //DEBUG_PRINT("backtrack", "Cumulative Relation: " << next_cell.pair_relation);
        //DEBUG_PRINT("backtrack", "Cumulative Query: " << next_cell.query_sequence);

        // Check if the we've reached the end of the subsequence, if not, add cell back to queue
        if(memo[next_cell.row_idx][next_cell.col_idx] != 0){backtrack_queue.push_back(next_cell);} 
        // Otherwise add the struct with our completed path to our results
        else {results.push_back(next_cell);}

    } // end while
    //cout << std::endl;
}

void BandedSmithWaterman::align(){
    init_matrix();
    #ifdef PRINT_MATRIX
        print_matrix();
    #endif
    score_matrix();
    backtrack();
    print_results();
}

void BandedSmithWaterman::print_results(){
    #ifdef PRINT_MATRIX
        cout << "[Scored Matrix]\n";
        print_matrix();
    #endif
    #ifdef PRINT_EXTRA
        cout << "[# Highest Scoring Sequences: " << results.size() << " Score: " << max_score << "]\n";
        cout << "[Sequence Pairing(s)]\n";
        cout << "====================\n";
    #else
        cout << max_score << "\n";
    #endif

    for(size_t idx = 0; idx < results.size(); idx++){
        #ifdef PRINT_EXTRA
            cout << "Sequence Pair: " << idx << "\n";
        #endif
        backtrack_info backtrack_record = results[idx];
        for(size_t base = 0; base < backtrack_record.query_sequence.size(); base++){
            cout << backtrack_record.query_sequence[base];
        }
        cout << "\n";
        for(size_t base = 0; base < backtrack_record.pair_relation.size(); base++){
            cout << backtrack_record.pair_relation[base];
        }
        cout << "\n";
        for(size_t base = 0; base < backtrack_record.reference_sequence.size(); base++){
            cout << backtrack_record.reference_sequence[base];
        }
        #ifdef PRINT_EXTRA
            cout << "\n====================\n";
        #else
            cout << "\n";
        #endif
    }

    //cout << std::endl;
}