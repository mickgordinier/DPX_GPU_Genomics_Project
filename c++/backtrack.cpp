#include  "backtrack.h"

void printMatrix(const int *memo, const int cols, const int rows){
    for(size_t row = 0; row < rows; row++){
        for(size_t col = 0; col < cols; col++){
            printf(" %4d ", memo[row*(cols) + col]);
        }
        printf("\n");
    }
}

void backtrackNW(const direction *backtrackMemo, const char *referenceString, const int referenceLength, const char *queryString,  const int queryLength){
    int currentMemoRow = queryLength;
    int currentMemoCol = referenceLength;

    std::string referenceSequence = "";
    std::string pairRelation = "";
    std::string querySequence = "";

    while ((currentMemoRow != 0) && (currentMemoCol != 0)) {
        
        // Determine the current cell's predecessor
        switch (backtrackMemo[(currentMemoRow * (queryLength+1)) + currentMemoCol]) {
            
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

    printf("%s\n", referenceSequence);
    printf("%s\n", pairRelation);
    printf("%s\n", querySequence);
    
}