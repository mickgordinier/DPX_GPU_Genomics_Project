#include "LinearSmithWaterman.h"

using std::cout;
using std::vector;
using std::max;
using std::string;

void LinearSmithWaterman::init_matrix(){
    cout << "[Initializing Matrix...]\n";
    size_t num_rows = query_str.size();
    size_t num_cols = reference_str.size();

    // Initialize scoring memo
    vector<int> temp_memo;
    temp_memo.resize(num_cols + 1, 0);
    memo.resize(num_rows + 1, temp_memo);

    // Initialize backtracking memo
    vector<direction> temp_backtrack;
    temp_backtrack.resize(num_cols, NONE);
    backtrack_memo.resize(num_rows, temp_backtrack);
}

void LinearSmithWaterman::print_matrix(){
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

void LinearSmithWaterman::score_matrix() {
    cout << "[Scoring...]\n";

    // Score the smith waterman matrix, keeping track of where a cell got its value
    for(size_t row_idx = 1; row_idx < memo.size(); row_idx++) {
        for(size_t col_idx = 1; col_idx < memo[0].size(); col_idx++) {
            int upper_score;
            int left_score;
            int corner_score;

            // Calculate potential gap scores
            upper_score = memo[row_idx - 1][col_idx] + gap_weight;
            left_score = memo[row_idx][col_idx - 1] + gap_weight;
            
            direction corner_direction;
            // If match, add the match score from the corner
            if (query_str[row_idx - 1] == reference_str[col_idx - 1]){
                corner_score = memo[row_idx - 1][col_idx - 1] + match_weight;
                corner_direction = CORNER_MATCH;
            // Otherwise, add the penalty score from the corner
            }
            else {
                corner_score = memo[row_idx - 1][col_idx - 1] + mismatch_weight;
                corner_direction = CORNER_MISMATCH;
            }

            // Find the largest of the 3 potential scores
            int temp_greatest = max(upper_score, max(left_score, corner_score));
            
            // Apply ReLU to score
            memo[row_idx][col_idx] = max(0, temp_greatest);
            
            // If the greatest amongst the 3 is still negative, no backtracking
            if (temp_greatest < 0) {continue;}
            // Otherwise, note where the cell needs to backtrack to
            else if(corner_score == memo[row_idx][col_idx]){backtrack_memo[row_idx-1][col_idx-1] = corner_direction;} 
            else if (left_score == memo[row_idx][col_idx]) {backtrack_memo[row_idx-1][col_idx-1] = LEFT_GAP;} 
            else {backtrack_memo[row_idx-1][col_idx-1] = UPPER_GAP;}// end if

        } // end for col_idx
    } // end for row_idx

    //cout << std::endl;
}

void LinearSmithWaterman::backtrack(){
    cout << "[Finding max cells...]\n";

    // We need to find the max score of the entire matrix
    // + note the cell(s) that have this value

    // Whenever we find a new max value we need to clear our backtracking queue
    // --> As a result, start from the bottom right so we hopefully 
    // end up clearing the queue fewer times
    for(size_t row_idx = memo.size() - 1; row_idx > 0; row_idx--) {
        for(size_t col_idx = memo[0].size() - 1; col_idx > 0; col_idx--) {
            int new_max_score = max(max_score, memo[row_idx][col_idx]);
            // Clear the queue if we have a new max val + store it
            if(new_max_score > max_score){
                backtrack_queue.clear();
                max_score = new_max_score;
                //DEBUG_PRINT("backtrack", "Clearing backtrack queue");
            }
            // Add the new starting cell to the backtracking queue
            if(memo[row_idx][col_idx] == max_score){
                backtrack_info new_info = {"", "", "", row_idx, col_idx};
                backtrack_queue.push_back(new_info);
                //DEBUG_PRINT("backtrack", "Adding [" << row_idx << ", " << col_idx << "] to backtrack queue with Val: " << max_score);
            }
        }
    }

    cout << "[Backtracking...]\n";
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

void LinearSmithWaterman::align(){
    init_matrix();
    #ifdef PRINT_MATRIX
        print_matrix();
    #endif
    score_matrix();
    backtrack();
    print_results();
}

void LinearSmithWaterman::print_results(){
    cout << "[Scored Matrix]\n";
    #ifdef PRINT_MATRIX
        print_matrix();
    #endif

    cout << "[# Highest Scoring Sequences: " << results.size() << " Score: " << max_score << "]\n";
    cout << "[Sequence Pairing(s)]\n";
    cout << "====================\n";

    for(size_t idx = 0; idx < results.size(); idx++){
        cout << "Sequence Pair: " << idx << "\n";
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
        cout << "\n====================\n";
    }

    //cout << std::endl;
}