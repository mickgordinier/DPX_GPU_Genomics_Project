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
    int tid = threadIdx.x;

    // Handling the matrix 0th row
    if (tid == 0) {
        for (int i = 0; i < referenceLength+1; ++i) {
            scoringMatrix[i] = i * gapWeight;
        }
    }

    int queryInsertionScore;
    int queryDeletionScore;
    int matchMismatchScore;
    int largestScore;

    int whileCount = 0;

    int rowIdx, matrixIdx, rowAboveIdx;
    char queryChar;

    while ((whileCount * BLOCK_SIZE) + 1) {

        rowIdx = (whileCount * BLOCK_SIZE) + tid + 1;

        matrixIdx = rowIdx * referenceLength;
        rowAboveIdx = ((rowIdx-1) * referenceLength) + 1;

        // Each thread initializing the 0th column of their row to 0
        if (rowIdx < queryLength) {
            scoringMatrix[matrixIdx++] = queryDeletionWeight * tid;
            queryChar = queryString[row_idx - 1]; 
        }

        // Need to ensure all threads in the block have written to their respective location before continuing
        __syncthreads();

        // Need to fill up all the rows
        // Starting on column idx 1
        for (int i = 1; i < referenceLength + BLOCK_SIZE - 1; ++i) {

            // On initialization, the lower threads need to wait for the upper thread to begin
            // if i < tid, we can assume the thread has not yet calculated the above value
            
            // At the end, we don't want the thread to do any more computation at the end of its row
            // Thus, we will have the thread stop updating the matrix once i >= (referenceLength + tid)
            if ((rowIdx < queryLength) && (i >= tid) && (i < (referenceLength + tid))) {

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

        ++whileCount;

    }
}


int main() {

    // TODO: NEED WAY TO RECEIVE USER INPUT
    // Need to be able to input reference and query strings
    // TEMPORARY HARD CODED STRING AND WEIGHT VALUES
    const string referenceString = "??????";
    const string queryString = "??????";

    const int matchWeight = 3;
    const int mismatchWeight = -2;
    const int gapWeight = -1;


    // Allocate device memory for matrices
    int *deviceScoringMatrix;
    direction *deviceBacktrackMatrix;
    cudaMalloc(&deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int));
    cudaMalloc(&deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction));


    // Need to launch kernel
    // Launching a kernel with 1 block with BLOCK_SIZE threads to populate scoring matrix
    needleman_wunsch_forward_pass_kernel<<<1, BLOCK_SIZE>>>(
        deviceScoringMatrix, deviceBacktrackMatrix,
        queryString, referenceString, queryString.size(), referenceString.size(), 
        matchWeight, mismatchWeight, gapWeight)


    // Allocate host memory for matrices
    // Allow for matrices to come from device -> host
    // Free up device memory
    int *hostScoringMatrix = new int[(referenceLength+1) * (queryLength+1)];
    int *hostBacktrackMatrix = new int[(referenceLength+1) * (queryLength+1)];

    cudaMemcpy(hostScoringMatrix, deviceScoringMatrix, (referenceLength+1) * (queryLength+1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hostBacktrackMatrix, deviceBacktrackMatrix, (referenceLength+1) * (queryLength+1) * sizeof(direction), cudaMemcpyDeviceToHost);

    cudaFree(deviceScoringMatrix);
    cudaFree(deviceBacktrackMatrix);


    // Perform backtracking on host and print results
    cout << "[Finding max cells...]\n";

    // We need to find the max score of the entire matrix
    // + note the cell(s) that have this value

    int currentMemoRow = hostBacktrackMatrix.size()-1;
    int currentMemoCol = hostBacktrackMatrix[0].size()-1;

    std::string referenceSequence = "";
    std::string pairRelation = "";
    std::string querySequence = "";

    while ((currentMemoRow != 0) && (currentMemoCol != 0)) {
        
        // Determine the current cell's predecessor
        switch (hostBacktrackMatrix[currentMemoRow][currentMemoCol]) {
            
            case MATCH:
                referenceSequence = referenceString[currentMemoCol-1] + referenceSequence;
                pairRelation = "*" + pairRelation;
                querySequence = queryString[currentMemoRow-1] + querySequence;
                --currentMemoRow;
                --currentMemoCol;
                break;
            // end if match

            case MISMATCH: 
                referenceSequence = referenceString[currentMemoCol-1] + referenceSequence;
                pairRelation = "|" + pairRelation;
                querySequence = queryString[currentMemoRow-1] + querySequence;
                --currentMemoRow;
                --currentMemoCol;
                break;
            // end if mismatch
            
            case QUERY_DELETION:
                referenceSequence = "_" + referenceSequence;
                pairRelation = " " + pairRelation;
                querySequence = queryString[currentMemoRow-1] + querySequence;
                --currentMemoRow;
                break;
            // end if query deletion
            
            case QUERY_INSERTION:
                referenceSequence = referenceString[currentMemoCol-1] + referenceSequence;
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

    cout << "[Needleman-Wunsch Score: " << hostScoringMatrix[hostScoringMatrix.size()-1][hostScoringMatrix[0].size()-1] << "]\n\n";
    cout << "[Sequence Pairing(s)]\n";
    cout << "\n====================\n";

    cout << referenceSequence << std::endl;
    cout << pairRelation << std::endl;
    cout << querySequence << std::endl;
    
    cout << "\n====================\n";
    cout << std::endl;

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