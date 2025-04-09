#include "LinearNeedlemanWunsch.h"
#include "FakeDPX.hpp"

using std::cout;
using std::vector;
using std::max;
using std::string;

// Enumerated predecessor directions for backtracking
enum direction {
    NONE,             // For initialization purposes
    MATCH,
    MISMATCH,
    QUERY_INSERTION,
    QUERY_DELETION
}; 

// Blocks are 1D with a size of the maximum 1024 threads
#define BLOCK_SIZE 1024

/*
    THINGS TO CONSIDER FOR OPTIMIZATION
    1. Complete removal of the scoring matrix altogether (Use of warp shuffling and shared memory)
    2. Using 16x2 DPX instructions to have a thread work on 2 cells concurrently


*/

// Device kernel that each thread will be executing to fill in its respective row
__global__ void
needleman_wunsch_forward_pass_kernel(int *scoringMatrix, direction *backtrackMatrix,
                        const string queryString, const string referenceString,
                        const int queryLength, const int referenceLength,
                        const int matchWeight, const int mismatchWeight, const int gapWeight)
{
    // Obtaining the 1D unique block and thread Id for the specific thread
    int bid = blockIdx.x;
    int tid = threadIdx.x;

    int queryInsertionScore;
    int queryDeletionScore;
    int matchMismatchScore;
    int largestScore;

    while () {

        // We will handle the 0th row in Device code
        int rowIdx = (whileCount * BLOCK_SIZE) + ((bid * BLOCK_SIZE) + tid) + 1;

        int matrixIdx = rowIdx * referenceLength;

        // Each thread initializing the 0th column of their row to 0
        scoringMatrix[matrixIdx++] = queryDeletionWeight * tid;

        int rowAboveIdx = ((rowIdx-1) * referenceLength) + 1;

        int queryInsertionScore;
        int queryDeletionScore;
        int matchMismatchScore;
        int largestScore;

        const char queryChar = queryString[row_idx - 1]; 


        // Need to fill up all the rows
        // Starting on column idx 1
        for (int i = 1; i < referenceLength + BLOCK_SIZE - 1; ++i) {

            // On initialization, the lower threads need to wait for the upper thread to begin
            // if i < tid, we can assume the thread has not yet calculated the above value
            
            // At the end, we don't want the thread to do any more computation at the end of its row
            // Thus, we will have the thread stop updating the matrix once i >= (referenceLength + tid)
            if ((i >= tid) & (i < (referenceLength + tid))) {

                corner_direction = NONE;

                const char referenceChar = referenceLength[i - tid];

                // If match, add the match score from the corner
                if (queryChar == referenceChar){
                    matchMismatchScore = scoringMatrix[rowAboveIdx-1] + matchWeight;
                    corner_direction = MATCH;
                }
                // Otherwise, add the penalty score from the corner
                else {
                    matchMismatchScore = scoringMatrix[rowAboveIdx-1] + mismatchWeight;
                    corner_direction = MISMATCH;
                }
                
                // Calculate potential gap scores
                queryDeletionScore = scoringMatrix[rowAboveIdx] + gapWeight;
                queryInsertionScore = scoringMatrix[matrixIdx-1] + gapWeight;
                
                largestScore = __vibmax_s32(queryDeletionScore, matchMismatchScore, &pred);
                if (pred) corner_direction = QUERY_DELETION;
                
                largestScore = __vibmax_s32(queryInsertionScore, largestScore, &pred);
                if (pred) corner_direction = QUERY_INSERTION;

                scoringMatrix[matrixIdx] = largestScore;
                backtrackMatrix[matrixIdx] = corner_direction;

                ++matrixIdx;
                ++rowAboveIdx;
            }

            // Need to ensure all threads in the block have written to their respective location before continuing
            __syncthreads();

        }
    }
}


int main() {

    // TODO: NEED WAY TO RECEIVE USER INPUT
    // Need to be able to input reference and query strings

    // Allocate device memory for matrices
    int *deviceScoringMatrix, deviceBacktrackMatrix;
    cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int));
    cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int));

    // Need to launch kernel

    // Need to perform bactracking





    // Allocate host memory for matrices
    // Allow for matrices to come from device -> host
    int *hostScoringMatrix = new int[(referenceLength+1) * (queryLength+1)];
    int *hostBacktrackMatrix = new int[(referenceLength+1) * (queryLength+1)];

    cudaMemcpy(hostScoringMatrix, deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hostBacktrackMatrix, deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(deviceScoringMatrix);
    cudaFree(deviceBacktrackMatrix);










    cout << "[Initializing Matrix...]\n";
    size_t num_rows = query_str.size();
    size_t num_cols = reference_str.size();

    // Initialize scoring memo
    // Need the additional 0th row and column for initial gap penalty scores
    vector<int> temp_memo;
    temp_memo.resize(num_cols + 1, 0);
    memo.resize(num_rows + 1, temp_memo);



}






void init_matrix(){
    cout << "[Initializing Matrix...]\n";
    size_t num_rows = query_str.size();
    size_t num_cols = reference_str.size();

    // Initialize scoring memo
    // Need the additional 0th row and column for initial gap penalty scores
    vector<int> temp_memo;
    temp_memo.resize(num_cols + 1, 0);
    memo.resize(num_rows + 1, temp_memo);

    // Initialize backtracking memo
    // As we are doing global alignment we need to be able to backtrack from 0th row and column
    vector<direction> temp_backtrack;
    temp_backtrack.resize(num_cols + 1, NONE);
    backtrack_memo.resize(num_rows + 1, temp_backtrack);

    // If query needed to do initial DELETIONS before entering reference
    // Can be thought of adding insertions to reference
    // As this is linear gap penalty, will be doing row_idx * gap_weight
    for (int row_idx = 1; row_idx < num_rows+1; ++row_idx) {
        memo[row_idx][0] = row_idx * gap_weight;
        backtrack_memo[row_idx][0] = QUERY_DELETION;
    }

    // If query needed to do initial INSERTIONS before entering reference
    // As this is linear gap penalty, will be doing row_idx * gap_weight
    for (int col_idx = 1; col_idx < num_cols+1; ++col_idx) {
        memo[0][col_idx] = col_idx * gap_weight;
        backtrack_memo[0][col_idx] = QUERY_INSERTION;
    }
}

void LinearNeedlemanWunsch::print_matrix(){
    
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

void LinearNeedlemanWunsch::score_matrix() {
    cout << "[Scoring...]\n";

    bool pred;

    // Score the smith waterman matrix, keeping track of where a cell got its value
    for(size_t row_idx = 1; row_idx < memo.size(); row_idx++) {
        for(size_t col_idx = 1; col_idx < memo[0].size(); col_idx++) {
            
            int queryInsertionScore;
            int queryDeletionScore;
            int matchMismatchScore;
            int largestScore;
            
            direction corner_direction = NONE;

            // If match, add the match score from the corner
            if (query_str[row_idx - 1] == reference_str[col_idx - 1]){
                matchMismatchScore = memo[row_idx - 1][col_idx - 1] + match_weight;
                corner_direction = MATCH;
            }
            // Otherwise, add the penalty score from the corner
            else {
                matchMismatchScore = memo[row_idx - 1][col_idx - 1] + mismatch_weight;
                corner_direction = MISMATCH;
            }
            
            // Calculate potential gap scores
            queryDeletionScore = memo[row_idx - 1][col_idx] + gap_weight;
            queryInsertionScore = memo[row_idx][col_idx - 1] + gap_weight;
            
            largestScore = FakeDPX::__vibmax_s32(queryDeletionScore, matchMismatchScore, &pred);
            if (pred) corner_direction = QUERY_DELETION;
            
            largestScore = FakeDPX::__vibmax_s32(queryInsertionScore, largestScore, &pred);
            if (pred) corner_direction = QUERY_INSERTION;

            memo[row_idx][col_idx] = largestScore;
            backtrack_memo[row_idx][col_idx] = corner_direction;

        } // end for col_idx
    } // end for row_idx

    cout << std::endl;
}

void LinearNeedlemanWunsch::backtrack(){
    cout << "[Finding max cells...]\n";

    // We need to find the max score of the entire matrix
    // + note the cell(s) that have this value

    int currentMemoRow = memo.size()-1;
    int currentMemoCol = memo[0].size()-1;

    std::string referenceSequence = "";
    std::string pairRelation = "";
    std::string querySequence = "";

    while ((currentMemoRow != 0) && (currentMemoCol != 0)) {
        
        // Determine the current cell's predecessor
        switch (backtrack_memo[currentMemoRow][currentMemoCol]) {
            
            case MATCH:
                referenceSequence = reference_str[currentMemoCol-1] + referenceSequence;
                pairRelation = "*" + pairRelation;
                querySequence = query_str[currentMemoRow-1] + querySequence;
                --currentMemoRow;
                --currentMemoCol;
                break;
            // end if match

            case MISMATCH: 
                referenceSequence = reference_str[currentMemoCol-1] + referenceSequence;
                pairRelation = "|" + pairRelation;
                querySequence = query_str[currentMemoRow-1] + querySequence;
                --currentMemoRow;
                --currentMemoCol;
                break;
            // end if mismatch
            
            case QUERY_DELETION:
                referenceSequence = "_" + referenceSequence;
                pairRelation = " " + pairRelation;
                querySequence = query_str[currentMemoRow-1] + querySequence;
                --currentMemoRow;
                break;
            // end if query deletion
            
            case QUERY_INSERTION:
                referenceSequence = reference_str[currentMemoCol-1] + referenceSequence;
                pairRelation = " " + pairRelation;
                querySequence = "_" + querySequence;
                --currentMemoCol;
                break;
            // end if query insertion
            
            default:
                exit(1);
            // end if upper gap

        } // end switch
    } // end while

    print_matrix();

    cout << "[Needleman-Wunsch Score: " << memo[memo.size()-1][memo[0].size()-1] << "]\n\n";
    cout << "[Sequence Pairing(s)]\n";
    cout << "\n====================\n";

    cout << referenceSequence << std::endl;
    cout << pairRelation << std::endl;
    cout << querySequence << std::endl;
    
    cout << "\n====================\n";
    cout << std::endl;
}

void LinearNeedlemanWunsch::align(){
    init_matrix();
    print_matrix();
    score_matrix();
    backtrack();
}

void LinearNeedlemanWunsch::print_results(){
//     cout << "[Scored Matrix]\n";
//     print_matrix();

//     cout << "[# Highest Scoring Sequences: " << results.size() << " Score: " << max_score << "]\n\n";
//     cout << "[Sequence Pairing(s)]\n";
//     cout << "\n====================\n";

//     for(size_t idx = 0; idx < results.size(); idx++){
//         cout << "Sequence Pair: " << idx << "\n";
//         backtrack_info backtrack_record = results[idx];
//         for(size_t base = 0; base < backtrack_record.query_sequence.size(); base++){
//             cout << backtrack_record.query_sequence[base];
//         }
//         cout << "\n";
//         for(size_t base = 0; base < backtrack_record.pair_relation.size(); base++){
//             cout << backtrack_record.pair_relation[base];
//         }
//         cout << "\n";
//         for(size_t base = 0; base < backtrack_record.reference_sequence.size(); base++){
//             cout << backtrack_record.reference_sequence[base];
//         }
//         cout << "\n====================\n";
//     }

//     cout << std::endl;
}


