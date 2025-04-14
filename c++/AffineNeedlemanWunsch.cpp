#include "AffineNeedlemanWunsch.h"
#include "FakeDPX.hpp"

using std::cout;
using std::vector;
using std::max;
using std::string;

// #define PRINT_MATRIX
// #define PRINT_EXTRA

void AffineNeedlemanWunsch::init_matrix(){
    
    #ifdef PRINT_EXTRA
        cout << "[Initializing Matrix...]\n";
    #endif

    size_t num_rows = query_str.size();
    size_t num_cols = reference_str.size();

    // Initialize scoring memo
    // Need the additional 0th row and column for initial gap penalty scores
    vector<int> tempRow;
    tempRow.resize(num_cols + 1, 0);
    scoringMemo.resize(num_rows + 1, tempRow);
    queryInsertionMemo.resize(num_rows + 1, tempRow);
    queryDeletionMemo.resize(num_rows + 1, tempRow);

    // Initialize backtracking memo
    // As we are doing global alignment we need to be able to backtrack from 0th row and column
    vector<directionMain> temp_backtrack;
    temp_backtrack.resize(num_cols + 1, NONE_MAIN);
    scoringBacktrack.resize(num_rows + 1, temp_backtrack);

    vector<directionIndel> temp_backtrack_2;
    temp_backtrack_2.resize(num_cols + 1, NONE_INDEL);
    queryInsertionBacktrack.resize(num_rows + 1, temp_backtrack_2);
    queryDeletionBacktrack.resize(num_rows + 1, temp_backtrack_2);

    // If query needed to do initial DELETIONS before entering reference
    // Can be thought of adding insertions to reference
    // As this is linear gap penalty, will be doing rowIdx * gap_weight
    for (int rowIdx = 1; rowIdx < num_rows+1; ++rowIdx) {
        scoringMemo[rowIdx][0] = gapOpenWeight + (rowIdx * gapExtendWeight);
        scoringBacktrack[rowIdx][0] = QUERY_DELETION;
    }

    // If query needed to do initial INSERTIONS before entering reference
    // As this is linear gap penalty, will be doing rowIdx * gap_weight
    for (int colIdx = 1; colIdx < num_cols+1; ++colIdx) {
        scoringMemo[0][colIdx] = gapOpenWeight + (colIdx * gapExtendWeight);
        scoringBacktrack[0][colIdx] = QUERY_INSERTION;
    }
}

void AffineNeedlemanWunsch::print_matrix(){
    
    size_t num_rows = query_str.size();
    size_t num_cols = reference_str.size();

    // Print parameters
    cout << "Reference: " << reference_str << " Size: " << num_cols << "\n";
    cout << "Query: " << query_str << " Size: " << num_rows << "\n";

    cout << "Matrix Dim: [ " << num_rows + 1 << " x " << num_cols + 1 << " ]\n\n";

    // Set the spacing for each element of the matrix
    int min_width = 2;

    cout << "SCORING MATRIX\n";

    // Print out reference str
    cout << "    " << std::setw(min_width) << " " << " ";
    for(size_t i = 0; i < reference_str.size(); i++){
        cout << "  " << std::setw(min_width) << reference_str[i] << " ";
    }
    cout << "\n";

    // Print rest of matrix
    for(size_t i = 0; i < scoringMemo.size(); i++){
        if(i == 0){
            cout << "  ";
        }
        else {
            cout << query_str[i-1] << " ";
        }
        cout << "[";
        for(size_t j = 0; j < scoringMemo[i].size(); j++){
            if(j != scoringMemo[i].size() - 1){
                cout << " " << std::setw(min_width) << scoringMemo[i][j] << ", ";
            }
            else {
                cout << " " << std::setw(min_width) << scoringMemo[i][j];
            }
        } //end for
        cout << "],\n";
    } //end for
    
    cout << "\n";

    cout << "QUERY INSERTION MATRIX\n";

    // Print out reference str
    cout << "    " << std::setw(min_width) << " " << " ";
    for(size_t i = 0; i < reference_str.size(); i++){
        cout << "  " << std::setw(min_width) << reference_str[i] << " ";
    }
    cout << "\n";

    // Print rest of matrix
    for(size_t i = 0; i < queryInsertionMemo.size(); i++){
        if(i == 0){
            cout << "  ";
        }
        else {
            cout << query_str[i-1] << " ";
        }
        cout << "[";
        for(size_t j = 0; j < queryInsertionMemo[i].size(); j++){
            if(j != queryInsertionMemo[i].size() - 1){
                cout << " " << std::setw(min_width) << queryInsertionMemo[i][j] << ", ";
            }
            else {
                cout << " " << std::setw(min_width) << queryInsertionMemo[i][j];
            }
        } //end for
        cout << "],\n";
    } //end for
    
    cout << "\n";

    cout << "QUERY DELETION MATRIX\n";

    // Print out reference str
    cout << "    " << std::setw(min_width) << " " << " ";
    for(size_t i = 0; i < reference_str.size(); i++){
        cout << "  " << std::setw(min_width) << reference_str[i] << " ";
    }
    cout << "\n";

    // Print rest of matrix
    for(size_t i = 0; i < queryDeletionMemo.size(); i++){
        if(i == 0){
            cout << "  ";
        }
        else {
            cout << query_str[i-1] << " ";
        }
        cout << "[";
        for(size_t j = 0; j < queryDeletionMemo[i].size(); j++){
            if(j != queryDeletionMemo[i].size() - 1){
                cout << " " << std::setw(min_width) << queryDeletionMemo[i][j] << ", ";
            }
            else {
                cout << " " << std::setw(min_width) << queryDeletionMemo[i][j];
            }
        } //end for
        cout << "],\n";
    } //end for

    cout << std::endl;

}

void AffineNeedlemanWunsch::score_matrix() {
    #ifdef PRINT_EXTRA
        cout << "[Scoring...]\n";
    #endif

    bool pred;

    // Score the smith waterman matrix, keeping track of where a cell got its value
    for(size_t rowIdx = 1; rowIdx < scoringMemo.size(); rowIdx++) {
        for(size_t colIdx = 1; colIdx < scoringMemo[0].size(); colIdx++) {
            
            int queryInsertionScore;
            int queryDeletionScore;
            int matchMismatchScore;
            int largestScore;

            // Handling scores of performing an query deletion at the end
            // Calculating best score of either creating or extending the deletion gap
            if (rowIdx == 1) {
                // PROBABLY CAN HANDLE ROW 1 DURING INITIALIZATION PHASE
                // Always assuming just opening new gap
                queryDeletionMemo[rowIdx][colIdx] = scoringMemo[rowIdx-1][colIdx] + gapOpenWeight + gapExtendWeight;
                queryDeletionBacktrack[rowIdx][colIdx] = GAP_OPEN;
            } else {
                queryDeletionMemo[rowIdx][colIdx] = FakeDPX::__vibmax_s32(
                    scoringMemo[rowIdx-1][colIdx] + gapOpenWeight + gapExtendWeight,  // Opening new gap at the end
                    queryDeletionMemo[rowIdx-1][colIdx] + gapExtendWeight,            // Extending current gap at end
                    &pred
                );
                queryDeletionBacktrack[rowIdx][colIdx] = pred ? GAP_OPEN : GAP_EXTEND;
            }

            // Handling scores of performing an query insertion at the end
            // Calculating best score of either creating or extending the insertion gap
            if (colIdx == 1) {
                // PROBABLY CAN HANDLE COL 1 DURING INITIALIZATION PHASE
                // Always assuming just opening new gap
                queryInsertionMemo[rowIdx][colIdx] = scoringMemo[rowIdx][colIdx-1] + gapOpenWeight + gapExtendWeight;
                queryInsertionBacktrack[rowIdx][colIdx] = GAP_OPEN;
            } else {
                queryInsertionMemo[rowIdx][colIdx] = FakeDPX::__vibmax_s32(
                    scoringMemo[rowIdx][colIdx-1] + gapOpenWeight + gapExtendWeight,  // Opening new gap at the end
                    queryInsertionMemo[rowIdx][colIdx-1] + gapExtendWeight,           // Extending current gap at end
                    &pred
                );
                queryInsertionBacktrack[rowIdx][colIdx] = pred ? GAP_OPEN : GAP_EXTEND;
            }

            // Calculating distance matrix element score
            directionMain corner_direction = NONE_MAIN;

            // If match, add the match score from the corner
            if (query_str[rowIdx - 1] == reference_str[colIdx - 1]){
                matchMismatchScore = scoringMemo[rowIdx - 1][colIdx - 1] + matchWeight;
                corner_direction = MATCH;
            }
            // Otherwise, add the penalty score from the corner
            else {
                matchMismatchScore = scoringMemo[rowIdx - 1][colIdx - 1] + mismatchWeight;
                corner_direction = MISMATCH;
            }
            
            largestScore = FakeDPX::__vibmax_s32(queryDeletionMemo[rowIdx][colIdx], matchMismatchScore, &pred);
            if (pred) corner_direction = QUERY_DELETION;
            
            largestScore = FakeDPX::__vibmax_s32(queryInsertionMemo[rowIdx][colIdx], largestScore, &pred);
            if (pred) corner_direction = QUERY_INSERTION;

            scoringMemo[rowIdx][colIdx] = largestScore;
            scoringBacktrack[rowIdx][colIdx] = corner_direction;

        } // end for colIdx
    } // end for rowIdx
}

void AffineNeedlemanWunsch::backtrack(){
    #ifdef PRINT_EXTRA
        cout << "[Finding max cells...]\n";
    #endif
    // We need to find the max score of the entire matrix
    // + note the cell(s) that have this value

    int currentMemoRow = scoringMemo.size()-1;
    int currentMemoCol = scoringMemo[0].size()-1;
    currentMatrixPosition currentMatrix = SCORING;

    std::string referenceSequence = "";
    std::string pairRelation = "";
    std::string querySequence = "";

    while ((currentMemoRow != 0) && (currentMemoCol != 0)) {

        if (currentMatrix == SCORING) {

            // Determine the current cell's predecessor
            switch (scoringBacktrack[currentMemoRow][currentMemoCol]) {
                
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
                    // referenceSequence = "_" + referenceSequence;
                    // pairRelation = " " + pairRelation;
                    // querySequence = query_str[currentMemoRow-1] + querySequence;
                    // --currentMemoRow;
                    currentMatrix = DELETION;
                    break;
                    // end if query deletion
                
                case QUERY_INSERTION:
                    // referenceSequence = reference_str[currentMemoCol-1] + referenceSequence;
                    // pairRelation = " " + pairRelation;
                    // querySequence = "_" + querySequence;
                    // --currentMemoCol;
                    currentMatrix = INSERTION;
                    break;
                    // end if query insertion
                
                default:
                    cout << "ERROR: MAIN BACKTRACK BAD\n";
                    exit(1);
                // end if upper gap

            } // end switch
        }  // end if in SCORING matrix

        else if (currentMatrix == INSERTION) {

            referenceSequence = reference_str[currentMemoCol-1] + referenceSequence;
            pairRelation = " " + pairRelation;
            querySequence = "_" + querySequence;
            --currentMemoCol;

            switch (queryInsertionBacktrack[currentMemoRow][currentMemoCol]) {
                
                // Meaning that we came from the scoring matrix
                case GAP_OPEN:
                    currentMatrix = SCORING;
                    break;

                case GAP_EXTEND:
                    currentMatrix = INSERTION;
                    break;

                default:
                    cout << "ERROR: INSERTION BACKTRACK BAD\n";
                    exit(1);

            } // end switch
        }  // end if in QUERY INSERTION matrix

        else if (currentMatrix == DELETION) {

            referenceSequence = "_" + referenceSequence;
            pairRelation = " " + pairRelation;
            querySequence = query_str[currentMemoRow-1] + querySequence;
            --currentMemoRow;

            switch (queryDeletionBacktrack[currentMemoRow][currentMemoCol]) {
                
                // Meaning that we came from the scoring matrix
                case GAP_OPEN:
                    currentMatrix = SCORING;
                    break;

                case GAP_EXTEND:
                    currentMatrix = DELETION;
                    break;

                default:
                    cout << "ERROR: DELETION BACKTRACK BAD\n";
                    exit(1);

            } // end switch
        }  // end if in QUERY DELETION matrix

        else {
            exit(1);
        }

    } // end while

    while (currentMemoRow > 0) {
        referenceSequence = "_" + referenceSequence;
        pairRelation = " " + pairRelation;
        querySequence = query_str[currentMemoRow-1] + querySequence;
        --currentMemoRow;
    }

    while (currentMemoCol > 0) {
        referenceSequence = reference_str[currentMemoCol-1] + referenceSequence;
        pairRelation = " " + pairRelation;
        querySequence = "_" + querySequence;
        --currentMemoCol;
    }

    #ifdef PRINT_EXTRA
        cout << "[Needleman-Wunsch Score: " << scoringMemo[scoringMemo.size()-1][scoringMemo[0].size()-1] << "]\n";
        cout << "[Sequence Pairing(s)]\n";
        cout << "====================\n";
    #else
        cout << scoringMemo[scoringMemo.size()-1][scoringMemo[0].size()-1] << "\n";
    #endif

    cout << referenceSequence << std::endl;
    cout << pairRelation << std::endl;
    cout << querySequence << std::endl;
    
    #ifdef PRINT_EXTRA
        cout << "====================\n";
        cout << std::endl;
    #endif
}

void AffineNeedlemanWunsch::align(){
    init_matrix();
    #ifdef PRINT_MATRIX
        print_matrix();
    #endif
    score_matrix();
    #ifdef PRINT_MATRIX
        print_matrix();
    #endif
    backtrack();
}

void AffineNeedlemanWunsch::print_results(){
}