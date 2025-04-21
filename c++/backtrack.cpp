#include  "backtrack.h"

void printMatrix(const int *memo, const int cols, const int rows){
    for(size_t row = 0; row < rows; row++){
        for(size_t col = 0; col < cols; col++){
            printf(" %4d ", memo[row*(cols) + col]);
        }
        printf("\n");
    }
}

void printBacktrackMatrix(const directionMain *memo, const int cols, const int rows){
    for(size_t row = 0; row < rows; row++){
        for(size_t col = 0; col < cols; col++){
            printf(" %4d ", memo[row*(cols) + col]);
        }
        printf("\n");
    }
}

void backtrackNW(const directionMain *backtrackMemo, const char *referenceString, const int referenceLength, const char *queryString,  const int queryLength){
    int numRows = queryLength + 1;
    int numCols = referenceLength + 1;
    int currentMemoRow = queryLength;
    int currentMemoCol = referenceLength;

    std::string referenceSequence = "";
    std::string pairRelation = "";
    std::string querySequence = "";

    while ((currentMemoRow != 0) || (currentMemoCol != 0)) {
        
        // Determine the current cell's predecessor
        switch (backtrackMemo[(currentMemoRow * numCols) + currentMemoCol]) {
            
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

    printf("%s\n", referenceSequence.c_str());
    printf("%s\n", pairRelation.c_str());
    printf("%s\n", querySequence.c_str());
    
}

void backtrackSW(const directionMain *backtrackMemo, const char *referenceString, const int referenceLength, const char *queryString,  const int queryLength, const int *scoreMatrix){
    int numRows = queryLength + 1;
    int numCols = referenceLength + 1;
    
    int currentMemoRow = 0;
    int currentMemoCol = 0;

    int maxScore = 0;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            if (scoreMatrix[(i * numCols) + j] > maxScore) {
                maxScore = scoreMatrix[(i * numCols) + j];
                currentMemoRow = i;
                currentMemoCol = j;
            }
        }
    }

    std::string referenceSequence = "";
    std::string pairRelation = "";
    std::string querySequence = "";

    while (currentMemoRow > 0 && currentMemoCol > 0 && scoreMatrix[(currentMemoRow * numCols) + currentMemoCol] != 0) {
        
        // Determine the current cell's predecessor
        switch (backtrackMemo[(currentMemoRow * numCols) + currentMemoCol]) {
            
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

    printf("%s\n", referenceSequence.c_str());
    printf("%s\n", pairRelation.c_str());
    printf("%s\n", querySequence.c_str());
    
}

void backtrackMultiNW(const directionMain *backtrackMemo, const char *referenceString, const int referenceLength, const char *queryString, const int queryLength, const int pairNum, const int score){
    int numRows = queryLength + 1;
    int numCols = referenceLength + 1;
    int currentMemoRow = queryLength;
    int currentMemoCol = referenceLength;

    std::string referenceSequence = "";
    std::string pairRelation = "";
    std::string querySequence = "";

    while ((currentMemoRow != 0) || (currentMemoCol != 0)) {
        
        // Determine the current cell's predecessor
        switch (backtrackMemo[(currentMemoRow * numCols) + currentMemoCol]) {
            
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
                printLock();
                printf("Exiting(1) backtrack: %d\n", pairNum);
                printUnlock();
                exit(1);
            // end if upper gap

        } // end switch
    } // end while

    printLock();
    printf("%d | %d\n", pairNum, score);
    printf("%s\n", referenceSequence.c_str());
    printf("%s\n", pairRelation.c_str());
    printf("%s\n", querySequence.c_str());
    printUnlock();
    
}

void backtrackANW(
    const directionMain *scoringBacktrack, 
    const directionIndel *queryInsertionBacktrack,
    const directionIndel *queryDeletionBacktrack,
    const char *referenceString, const int referenceLength, 
    const char *queryString,  const int queryLength)
{
    int numRows = queryLength + 1;
    int numCols = referenceLength + 1;

    int currentMemoRow = queryLength;
    int currentMemoCol = referenceLength;

    currentMatrixPosition currentMatrix = SCORING;

    std::string referenceSequence = "";
    std::string pairRelation = "";
    std::string querySequence = "";

    while ((currentMemoRow != 0) && (currentMemoCol != 0)) {

        if (currentMatrix == SCORING) {

            // Determine the current cell's predecessor
            switch (scoringBacktrack[(currentMemoRow * numCols) + currentMemoCol]) {
                
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
                    // referenceSequence = "_" + referenceSequence;
                    // pairRelation = " " + pairRelation;
                    // querySequence = queryString[currentMemoRow-1] + querySequence;
                    // --currentMemoRow;
                    currentMatrix = DELETION;
                    break;
                    // end if query deletion
                
                case QUERY_INSERTION:
                    // referenceSequence = referenceString[currentMemoCol-1] + referenceSequence;
                    // pairRelation = " " + pairRelation;
                    // querySequence = "_" + querySequence;
                    // --currentMemoCol;
                    currentMatrix = INSERTION;
                    break;
                    // end if query insertion
                
                default:
                    exit(1);
                // end if upper gap

            } // end switch
        }  // end if in SCORING matrix

        else if (currentMatrix == INSERTION) {

            switch (queryInsertionBacktrack[(currentMemoRow * numCols) + currentMemoCol]) {
                
                // Meaning that we came from the scoring matrix
                case GAP_OPEN:
                    currentMatrix = SCORING;
                    break;

                case GAP_EXTEND:
                    currentMatrix = INSERTION;
                    break;

                default:
                    exit(1);

            } // end switch

            referenceSequence = referenceString[currentMemoCol-1] + referenceSequence;
            pairRelation = " " + pairRelation;
            querySequence = "_" + querySequence;
            --currentMemoCol;

        }  // end if in QUERY INSERTION matrix

        else if (currentMatrix == DELETION) {

            switch (queryDeletionBacktrack[(currentMemoRow * numCols) + currentMemoCol]) {
                
                // Meaning that we came from the scoring matrix
                case GAP_OPEN:
                    currentMatrix = SCORING;
                    break;

                case GAP_EXTEND:
                    currentMatrix = DELETION;
                    break;

                default:
                    exit(1);

            } // end switch

            referenceSequence = "_" + referenceSequence;
            pairRelation = " " + pairRelation;
            querySequence = queryString[currentMemoRow-1] + querySequence;
            --currentMemoRow;
            
        }  // end if in QUERY DELETION matrix

        else {
            exit(1);
        }

    } // end while

    while (currentMemoRow > 0) {
        referenceSequence = "_" + referenceSequence;
        pairRelation = " " + pairRelation;
        querySequence = queryString[currentMemoRow-1] + querySequence;
        --currentMemoRow;
    }

    while (currentMemoCol > 0) {
        referenceSequence = referenceString[currentMemoCol-1] + referenceSequence;
        pairRelation = " " + pairRelation;
        querySequence = "_" + querySequence;
        --currentMemoCol;
    }

    printf("%s\n", referenceSequence.c_str());
    printf("%s\n", pairRelation.c_str());
    printf("%s\n", querySequence.c_str());
}