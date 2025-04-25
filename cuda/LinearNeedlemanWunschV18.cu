#include <stdio.h>  // For printf()
#include <cstring> // Determining length of string
#include "../c++/parseInput.h"
#include "../c++/backtrack.h"
#include "../c++/timing.h"

// Blocks are 1D with a size of the 32 threads (For 1 warp)
#define BLOCK_SIZE 32
#define BATCH_SIZE 10

struct s16x2 {
    int16_t A;
    int16_t B;
};

__device__ uint32_t pack_s16x2(uint16_t valA, uint16_t valB){
    return ((static_cast<uint32_t>(valA) << 16) | static_cast<uint32_t>(valB));
}

__device__ __host__ s16x2 unpack_s16x2(uint32_t val){
    return {static_cast<int16_t>(val >> 16), static_cast<int16_t>(val & 0xFFFF)};
}

__device__ void backtracking(char *backtrackStringsRet, directionMain *backtrackMatrix, const char * referenceString, const char * queryString, int *stringSpacing, int numRows, int numCols, int referenceStrIdx, int alignmentStrIdx, int queryStrIdx){
    backtrackStringsRet[referenceStrIdx] = '\0';
    backtrackStringsRet[alignmentStrIdx] = '\0';
    backtrackStringsRet[queryStrIdx] = '\0';

    int currentMemoRow = numRows - 1;
    int currentMemoCol = numCols - 1;

    while ((currentMemoRow != 0) && (currentMemoCol != 0)) {

        referenceStrIdx--;
        alignmentStrIdx--;
        queryStrIdx--;

        int newBackMatrixCol = (currentMemoRow-1) % BLOCK_SIZE;  // Also technically tid (loosely)
        int stripeIdx = (currentMemoRow-1) / BLOCK_SIZE;
        int newBackMatrixRowOffset = (currentMemoCol-1) + newBackMatrixCol;
        int newBackMatrixRow = (stripeIdx * (numCols + BLOCK_SIZE - 1)) + newBackMatrixRowOffset;

        int matrixIdx = (newBackMatrixRow * BLOCK_SIZE) + newBackMatrixCol;
        
        // Determine the current cell's predecessor
        switch (backtrackMatrix[matrixIdx]) {
            
            case MATCH:
                backtrackStringsRet[referenceStrIdx] = referenceString[currentMemoCol-1];
                backtrackStringsRet[alignmentStrIdx] = '*';
                backtrackStringsRet[queryStrIdx] = queryString[currentMemoRow-1];
                --currentMemoRow;
                --currentMemoCol;
                break;
            // end if match

            case MISMATCH: 
                backtrackStringsRet[referenceStrIdx] = referenceString[currentMemoCol-1];
                backtrackStringsRet[alignmentStrIdx] = '|';
                backtrackStringsRet[queryStrIdx] = queryString[currentMemoRow-1];
                --currentMemoRow;
                --currentMemoCol;
                break;
            // end if mismatch
            
            case QUERY_DELETION:
                backtrackStringsRet[referenceStrIdx] = '_';
                backtrackStringsRet[alignmentStrIdx] = ' ';
                backtrackStringsRet[queryStrIdx] = queryString[currentMemoRow-1];
                --currentMemoRow;
                break;
            // end if query deletion
            
            case QUERY_INSERTION:
                backtrackStringsRet[referenceStrIdx] = referenceString[currentMemoCol-1];
                backtrackStringsRet[alignmentStrIdx] = ' ';
                backtrackStringsRet[queryStrIdx] = '_';
                --currentMemoCol;
                break;
            // end if query insertion
            
            default:
                printf("ERROR: numRows: %d, numCols %d, currentRow: %d, currentCol: %d Element Idx %d\n", numRows, numCols, currentMemoRow, currentMemoCol, matrixIdx);
                return;
            // end if upper gap

        } // end switch
    } // end while

    while (currentMemoRow != 0) {
        referenceStrIdx--;
        alignmentStrIdx--;
        queryStrIdx--;
        backtrackStringsRet[referenceStrIdx] = '_';
        backtrackStringsRet[alignmentStrIdx] = ' ';
        backtrackStringsRet[queryStrIdx] = queryString[currentMemoRow-1];
        --currentMemoRow;
    }

    while (currentMemoCol != 0) {
        referenceStrIdx--;
        alignmentStrIdx--;
        queryStrIdx--;
        backtrackStringsRet[referenceStrIdx] = referenceString[currentMemoCol-1];
        backtrackStringsRet[alignmentStrIdx] = ' ';
        backtrackStringsRet[queryStrIdx] = '_';
        --currentMemoCol;
    }

    stringSpacing[blockIdx.x] = referenceStrIdx;
}

__global__ void 
needleman_wunsch_kernel(
    uint32_t *similarityScores,
    int *stringSpacing,
    directionMain *batchBacktrackMatrices,
    int *batchIndices,
    char *backtrackStringsRet, 
    const char *allSequences, const seqPair *allSequenceInfo,
    const int16_t matchWeight, const int16_t mismatchWeight, const int16_t gapWeight,
    const int startingSequenceIdx, const int stringLengthMax)
{

    const int tid = threadIdx.x;

    extern __shared__ uint32_t warpEdgeScore[]; 

    // We are launching multiple blocks, each of a warp of threads
    // Each block handles their own sequence alignment
    // We index into the array to obtain the strings and length
    
    /* each block will compute 2 alignments */
    int batchIndexA = batchIndices[blockIdx.x*2];
    int batchIndexB = batchIndices[blockIdx.x*2 + 1];
    directionMain *backtrackMatrixA = batchBacktrackMatrices + batchIndexA;
    directionMain *backtrackMatrixB = batchBacktrackMatrices + batchIndexB;

    const int sequenceIdxA = startingSequenceIdx + 2*blockIdx.x;
    const seqPair sequenceInfoA = allSequenceInfo[sequenceIdxA];
    const int queryLengthA = sequenceInfoA.querySize;
    const int referenceLengthA = sequenceInfoA.referenceSize;

    const int sequenceIdxB = startingSequenceIdx + 2*blockIdx.x + 1;
    const seqPair sequenceInfoB = allSequenceInfo[sequenceIdxB];
    const int queryLengthB = sequenceInfoB.querySize;
    const int referenceLengthB = sequenceInfoB.referenceSize;
    
    const char *queryStringA = allSequences + sequenceInfoA.queryIdx;
    const char *referenceStringA = allSequences + sequenceInfoA.referenceIdx;

    const char *queryStringB = allSequences + sequenceInfoB.queryIdx;
    const char *referenceStringB = allSequences + sequenceInfoB.referenceIdx;

    // The matrices are of size (queryLength + 1) * (referenceLength + 1)
    const int numRowsA = sequenceInfoA.querySize + 1;
    const int numColsA = sequenceInfoA.referenceSize + 1;

    const int numRowsB = sequenceInfoB.querySize + 1;
    const int numColsB = sequenceInfoB.referenceSize + 1;

    int16_t finalScoreA = 0, finalScoreB = 0;

    if (tid == 0 && sequenceIdxA == 0){
    }
    
    /* use the longer of the two alignments to determine loop bounds */
    const int numRows = max(numRowsA, numRowsB);
    const int numCols = max(numColsA, numColsB);
    
    if (tid == 0 && sequenceIdxA == 0){
        printf("sequenceIdx %d numRowsA %d numColsA %d %s %s\n", sequenceIdxA, numRowsA, numColsA, queryStringA, referenceStringA);
        printf("sequenceIdx %d numRowsB %d numColsB %d\n", sequenceIdxB, numRowsB, numColsB);
        printf("maxRows %d maxCols %d\n", numRows, numCols);
    }

    /* --- (BEGIN) COMPUTING SCORES -- */
    uint32_t leftDiag = pack_s16x2(gapWeight*tid, gapWeight*tid);
    uint32_t left = pack_s16x2(gapWeight*(tid+1), gapWeight*(tid+1));
    uint32_t up = pack_s16x2(gapWeight*(tid+1), gapWeight*(tid+1));

    // if (tid == 0){
    //     printf("left unpacked vals: %x, packed: %x\n", gapWeight*(tid+1), left);
    //     printf("leftDiag unpacked vals: %x, packed: %x\n", gapWeight*tid, leftDiag);
    //     printf("up unpacked vals: %x, packed: %x\n", gapWeight*(tid+1), up);
    // }

    int stripeStartIdx = 0;

    char queryCharA = '\0', queryCharB = '\0';
    char referenceCharA = '\0', referenceCharB = '\0';

    // if (tid == 0){
    //     printf("seqIdxA %d | queryA %s | referenceA %s | numRowsA %d | numColsA %d | seqIdxB %d |queryB %s | referenceB %s | numRowsB %d | numColsB %d\n",
    //         sequenceIdxA, queryStringA, referenceStringA, numRowsA, numColsA, sequenceIdxB, queryStringB, referenceStringB, numRowsB, numColsB);
    // }

    // Going through all of the rows each thread has to do
    for (int stripeStart = 1; stripeStart < numRows; stripeStart+=BLOCK_SIZE){

        int row = stripeStart + tid;
        int32_t largestScore;

        /* threads outside of bounds should abort */
        if (row >= numRows) return;

        leftDiag = pack_s16x2(gapWeight*(row - 1), gapWeight*(row - 1));
        left = pack_s16x2(gapWeight*row, gapWeight*row);
        
        if ((row - 1) < queryLengthA){
            queryCharA = queryStringA[row-1];
        } else {
            queryCharA = '\0';
        }
        
        if ((row - 1) < queryLengthB){
            queryCharB = queryStringB[row - 1];
        } else {
            queryCharB = '\0';
        }
        
        // if (sequenceIdxA == 0) {
        //     printf("tid: %d, stripeStart: %d,  row: %d, leftDiag: %X, left %X, queryCharA: %c, queryCharB: %c\n", tid, stripeStart, row, leftDiag, left, queryCharA, queryCharB);
        // }

        for (int col = 1; col < (numCols+BLOCK_SIZE); ++col){
            int adj_col = col - tid;

            if (row == 1){
                leftDiag = pack_s16x2(gapWeight*(adj_col - 1), gapWeight*(adj_col - 1));
                up = pack_s16x2(gapWeight*adj_col, gapWeight*adj_col);
            }

            /* for all but the first stripe, t0 must grab its diagonal and upper values from t31 */
            if (stripeStart > 1 && tid == 0 && adj_col < numCols){
                up = warpEdgeScore[adj_col];
                leftDiag = (adj_col == 1) ? pack_s16x2(gapWeight*(row - 1), gapWeight*(row - 1)) : warpEdgeScore[adj_col - 1];
            }

            
            if (adj_col > 0 && adj_col < numCols){
                largestScore = pack_s16x2(0, 0);
                
                if (adj_col-1 < referenceLengthA){
                    referenceCharA = referenceStringA[adj_col-1];
                }
                
                if (adj_col-1 < referenceLengthB){
                    referenceCharB = referenceStringB[adj_col-1];
                }
                
                // if (tid == 0){
                //     printf("sequenceIdx %d | row %d adj_col %d | left %x | leftDiag %x | up %x\n", sequenceIdxA, row, adj_col, left, leftDiag, up);
                // }
                
                directionMain cornerDirectionA = NONE_MAIN;
                bool isMatchA = (queryCharA == referenceCharA);
                cornerDirectionA = isMatchA ? MATCH : MISMATCH;
                
                directionMain cornerDirectionB = NONE_MAIN;
                bool isMatchB = (queryCharB == referenceCharB);
                cornerDirectionB = isMatchB ? MATCH : MISMATCH;
                
                /* compute score for first alignment */
                s16x2 leftDiagUnpacked = unpack_s16x2(leftDiag);
                s16x2 leftUnpacked = unpack_s16x2(left);
                s16x2 upUnpacked = unpack_s16x2(up);

                // if (sequenceIdxA == 0 && row == 4 && adj_col == 2) {
                //     printf("Is A Match: %d\n", isMatchA);
                // }
                
                int16_t matchMismatchScoreA = isMatchA ? leftDiagUnpacked.A + matchWeight : leftDiagUnpacked.A + mismatchWeight;
                int16_t queryDeletionScoreA = upUnpacked.A + gapWeight;
                int16_t queryInsertionScoreA = leftUnpacked.A + gapWeight;
                
                /* compute score for second alignment */
                int16_t matchMismatchScoreB = isMatchB ? leftDiagUnpacked.B + matchWeight : leftDiagUnpacked.B + mismatchWeight;
                int16_t queryDeletionScoreB = upUnpacked.B + gapWeight;
                int16_t queryInsertionScoreB = leftUnpacked.B + gapWeight;
                
                uint32_t queryDeletionScore = pack_s16x2(queryDeletionScoreA, queryDeletionScoreB);
                uint32_t matchMismatchScore = pack_s16x2(matchMismatchScoreA, matchMismatchScoreB);
                uint32_t queryInsertionScore = pack_s16x2(queryInsertionScoreA, queryInsertionScoreB);
                    
                bool predA, predB;
                largestScore = __vibmax_s16x2(queryDeletionScore, matchMismatchScore, &predA, &predB);
                if (predA) cornerDirectionA = QUERY_DELETION;
                if (predB) cornerDirectionB = QUERY_DELETION;
                
                largestScore = __vibmax_s16x2(queryInsertionScore, largestScore, &predA, &predB);
                if (predA) cornerDirectionA = QUERY_INSERTION;
                if (predB) cornerDirectionB = QUERY_INSERTION;
                    
                /* 
                backtracking matrix is now padded for coalesced global memory accesses, so that
                all threads now write to elements in the same row, eg:
                t0  -   -   -   -   ... -   -
                t0  t1  -   -   -   ... -   -
                t0  t1  t2  -   -   ... -   - 
                ...
                -   -   -   -   -   ... t30 t31 
                -   -   -   -   -   ... -   t31
                */
                
                /*
                matrix A: 5x7 
                matrix B: 3x7
                
                STRIPE 0
                
                col = 3, tid = 2
                idxA: 0 + (col - 1)*4 + tid = 2*4 + 2 = 10
                idxB: 2*4 + 2 = 10
                
                
                stripe 1
                col = 6, tid = 0
                idxA: 1*(4 + 7 - 1)*4 + (5 * 4) + 0 = 60
                
                
                */
                int matrixIdxA = (stripeStartIdx * (BLOCK_SIZE+numColsA-1) * BLOCK_SIZE) + ((col - 1) * BLOCK_SIZE) + tid;
                int matrixIdxB = (stripeStartIdx * (BLOCK_SIZE+numColsB-1) * BLOCK_SIZE) + ((col - 1) * BLOCK_SIZE) + tid;
                
                
                if (adj_col < numColsA && row < numRowsA){
                    backtrackMatrixA[matrixIdxA] = cornerDirectionA;
                }
                
                if (adj_col < numColsB && row < numRowsB){
                    backtrackMatrixB[matrixIdxB] = cornerDirectionB;
                }
                
                // if (sequenceIdxA == 0 && adj_col == 188) {
                //     printf("(%c %c) tid: %d, row: %d, adj_col: %d, (LD, U, L) = (%x, %x, %x), Score = %x\n", queryCharA, referenceCharA, tid, row, adj_col, leftDiag, up, left, largestScore);
                // }

                left = largestScore;
                
                /* last thread in warp stores its scores in shared memory for t0 to access */
                if (tid == (BLOCK_SIZE - 1)){
                    warpEdgeScore[adj_col] = largestScore;
                }
                

                leftDiag = up;
            }
            

            /*  top value for thread n + 1 is thread n's largestScore (just calculated value)*/
            up = __shfl_up_sync(0xffffffff, largestScore, 1);

            // if (sequenceIdxA == 0){
            //     s16x2 largestScoreUnpacked = unpack_s16x2(largestScore);
            //     printf("tid %d row %d adj_col %d largestScore %x largestScoreA %d largestScoreB %d\n", tid, row, adj_col, largestScore, largestScoreUnpacked.A, largestScoreUnpacked.B);
            // }
        }

        if (row == numRowsA-1) {
            /* final score for alignment A is the upper 16 bits of the largestScore */
            finalScoreA = static_cast<int16_t>((largestScore >> 16));
            if (sequenceIdxA == 0){
                printf("final score A: %d\n", finalScoreA);
            } 
        }

        if (row == numRowsB-1){
            /* final score for alignment B is the lower 16 bits of the largestScore */
            finalScoreB = static_cast<int16_t>(largestScore & 0xFFFF);
            if (sequenceIdxA == 0){
                printf("final score B: %d\n", finalScoreB);
            } 
        }

        ++stripeStartIdx;
    }

    if (tid == 0) {
        /* NB: similarityScores contains uint32_t -- depending on host to split scores properly between two alignments */
        similarityScores[blockIdx.x] = pack_s16x2(finalScoreA, finalScoreB);
    }

    /* --- (END) POPULATING THE SCORING MATRIX -- */

    /* --- (BEGIN) DETERMINING BACKTRACKING -- */

    // Starting at the end
    if (tid == 0) {
        int referenceStrIdx = (stringLengthMax * 3) * 2 * blockIdx.x + (stringLengthMax-1);
        int alignmentStrIdx = referenceStrIdx + stringLengthMax;
        int queryStrIdx = alignmentStrIdx + stringLengthMax;

        backtracking(backtrackStringsRet, backtrackMatrixA, referenceStringA, queryStringA, stringSpacing, numRowsA, numColsA, referenceStrIdx, alignmentStrIdx, queryStrIdx);

        referenceStrIdx = (stringLengthMax * 3) * (2 * blockIdx.x + 1) + (stringLengthMax-1);
        alignmentStrIdx = referenceStrIdx + stringLengthMax;
        queryStrIdx = alignmentStrIdx + stringLengthMax;

        backtracking(backtrackStringsRet, backtrackMatrixB, referenceStringB, queryStringB, stringSpacing, numRowsB, numColsB, referenceStrIdx, alignmentStrIdx, queryStrIdx);
    }
}


void
handleErrs(
    cudaError_t err,
    const char *errMsg) 
{
    if (err != cudaSuccess) {
        printf(errMsg);
        printf("CUDA ERROR: %s\n", cudaGetErrorString(err));
        exit(1);
    }
}


int main(int argc, char *argv[]) {

    // Print some cuda details
    printf("[Cuda Details]\n");
    int deviceCount;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess) {
        printf("FAILED TO GET DEVICE COUNT\n");
        printf("CUDA test kernel error: %s\n", cudaGetErrorString(err));
        exit(1);
    }

    printf("Device count: %d\n", deviceCount);
    int device = 0;
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    printf("Device %d has compute capability %d.%d.\n",
           device, deviceProp.major, deviceProp.minor);
    printf("Concurrent kernels?: %d\n\n", deviceProp.concurrentKernels);

    // Check that YOU use it correctly
    if (argc < 2) {
		fprintf(stderr, "usage: main -pairs <InSeqFile> -match <matchWeight> -mismatch <mismatchWeight> -gap <gapWeight> \n");
		exit(EXIT_FAILURE);
    }
	
    // Get args
    char *pairFileName;
    int16_t matchWeight     = 3;
    int16_t mismatchWeight  = -1;
    int16_t gapWeight       = -2;
    if(strcmp(argv[1], "-pairs") == 0) {
        pairFileName = argv[2];
    }
    if(argc > 3 && strcmp(argv[3], "-match") == 0) {
        matchWeight = atoi(argv[4]);
    }
    if(argc > 5 && strcmp(argv[5], "-mismatch") == 0) {
        mismatchWeight = atoi(argv[6]);
    }
    if(argc > 7 && strcmp(argv[7], "-gap") == 0) {
        gapWeight = atoi(argv[8]);
    }

    // Parse input file
    printf("Parsing input file: %s\n", pairFileName);
    inputInfo fileInfo;
    seqPair* allSequenceInfo;
    char* sequences;
    fileInfo = parseInput(pairFileName, allSequenceInfo, sequences);
    printf("Num Pairs: %d\n\n", fileInfo.numPairs);

    // Start timer
    uint64_t kernel_time = 0;
    uint64_t memalloc_time = 0;
    uint64_t backtracking_time = 0;
    uint64_t printing_time = 0;
    uint64_t start_time = start_timer();
        
    // Copy over the sequences
    char* deviceSequences;
    seqPair *deviceAllSequenceInfo;

    uint64_t start_memalloc = get_time();
    handleErrs(
        cudaMalloc(&deviceSequences, (fileInfo.numBytes) * sizeof(char)),
        "FAILED TO ALLOCATE MEMORY FOR ALL SEQUENCES\n"
    );

    handleErrs(
        cudaMemcpy(deviceSequences, sequences, (fileInfo.numBytes) * sizeof(char), cudaMemcpyHostToDevice),
        "FAILED TO COPY MEMORY FOR ALL SEQUENCES\n"
    );

    handleErrs(
        cudaMalloc(&deviceAllSequenceInfo, (fileInfo.numPairs) * sizeof(seqPair)),
        "FAILED TO ALLOCATE MEMORY FOR ALL SEQUENCES\n"
    );

    handleErrs(
        cudaMemcpy(deviceAllSequenceInfo, allSequenceInfo, (fileInfo.numPairs) * sizeof(seqPair), cudaMemcpyHostToDevice),
        "FAILED TO COPY MEMORY FOR ALL SEQUENCES\n"
    );

    /* 
    store all backtracking matrices for a batch in one matrix - each warp will index in via index array
    */
    int *deviceBacktrackingIndices;
    int *hostBacktrackingIndices = (int *)malloc(BATCH_SIZE * sizeof(int));
    
    handleErrs(
        cudaMalloc(&deviceBacktrackingIndices, BATCH_SIZE * sizeof(int)),
        "FAILED TO ALLOCATE MEMORY TO deviceBacktrackingIndices\n"
    );
    
    uint32_t *deviceSimilarityScores;
    uint32_t *hostSimilarityScores = (uint32_t*)malloc(BATCH_SIZE * sizeof(uint32_t));

    handleErrs(
        cudaMalloc(&deviceSimilarityScores, BATCH_SIZE * sizeof(uint32_t)),
        "FAILED TO ALLOCATE MEMORY TO deviceSimilarityScores\n"
    );
    
    int *deviceStringSpacing;
    int *hostStringSpacing = (int*)malloc(BATCH_SIZE * sizeof(int));

    handleErrs(
        cudaMalloc(&deviceStringSpacing, BATCH_SIZE * sizeof(int)),
        "FAILED TO ALLOCATE MEMORY TO deviceSimilarityScores\n"
    );

    memalloc_time += get_time() - start_memalloc;

    // Run the kernel on every sequence
    for(size_t sequenceIdx = 0; sequenceIdx < fileInfo.numPairs; sequenceIdx+=BATCH_SIZE){
        start_memalloc = get_time();

        int largestReferenceLength = 0;
        int largestQueryLength = 0;

        /* first warp's starting index is 0 */
        uint64_t batchMatrixSize = 0;
        hostBacktrackingIndices[0] = 0;
        for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; ++i) {
            const int queryLength = allSequenceInfo[i].querySize;
            const int referenceLength = allSequenceInfo[i].referenceSize;

            largestReferenceLength = max(largestReferenceLength, referenceLength);
            largestQueryLength = max(largestQueryLength, queryLength);

            int numFullStripes = queryLength / BLOCK_SIZE;
            bool isQueryLeftover = (queryLength % BLOCK_SIZE != 0);

            int rowsPerStripe =  referenceLength + BLOCK_SIZE - 1;
            
            int totalRowsForNewBacktrack = (numFullStripes + isQueryLeftover) * rowsPerStripe;

            /* make sure we don't go over the end of the array */
            batchMatrixSize += (totalRowsForNewBacktrack * BLOCK_SIZE);
            if ((i - sequenceIdx) < (BATCH_SIZE - 1)){
                hostBacktrackingIndices[i-sequenceIdx + 1] = batchMatrixSize;
            }
        }

        /* copy backtracking indices to device */
        handleErrs(
            cudaMemcpy(deviceBacktrackingIndices, hostBacktrackingIndices, BATCH_SIZE * sizeof(int), cudaMemcpyHostToDevice),
            "FAILED TO COPY MEMORY FOR deviceBacktrackingIndices\n"
        );

        /* allocate device mem for all backtracking matrices */
        directionMain *deviceMatricesAll;
        handleErrs(
            cudaMalloc(&deviceMatricesAll, batchMatrixSize*sizeof(directionMain)),
            "FAILED TO ALLOCATE MEMORY TO deviceMatricesAll\n"
        );

        int stringLengthMax = (largestReferenceLength+largestQueryLength+1);

        char *deviceBacktrackStringRet;
    
        handleErrs(
            cudaMalloc(&deviceBacktrackStringRet, (stringLengthMax * 3) * BATCH_SIZE * sizeof(char)),
            "FAILED TO ALLOCATE MEMORY TO BACKTRACKING STRINGS\n"
        );

        memalloc_time += get_time() - start_memalloc;

        uint64_t start_kernel = get_time();
        // Need to launch kernel
        int smem_size = (largestReferenceLength + 1) * sizeof(uint32_t);
        needleman_wunsch_kernel<<<BATCH_SIZE, BLOCK_SIZE, smem_size>>>(
            deviceSimilarityScores,
            deviceStringSpacing,
            deviceMatricesAll,
            deviceBacktrackingIndices,
            deviceBacktrackStringRet,
            deviceSequences, deviceAllSequenceInfo,
            matchWeight, mismatchWeight, gapWeight,
            sequenceIdx, stringLengthMax
        );
        
        // Wait for kernel to finish
        handleErrs(
            cudaDeviceSynchronize(),
            "SYNCHRONIZATION FAILED\n"
        );
        kernel_time += get_time() - start_kernel;

        start_memalloc = get_time();

        handleErrs(
            cudaMemcpy(hostSimilarityScores, deviceSimilarityScores, BATCH_SIZE * sizeof(uint32_t), cudaMemcpyDeviceToHost),
            "FAILED TO COPY SIMILARITY SCORES FROM DEVICE --> HOST\n"
        );

        handleErrs(
            cudaMemcpy(hostStringSpacing, deviceStringSpacing, BATCH_SIZE * sizeof(int), cudaMemcpyDeviceToHost),
            "FAILED TO COPY SIMILARITY SCORES FROM DEVICE --> HOST\n"
        );

        char *hostBacktrackingStringRet = (char *)malloc(stringLengthMax * 3 * BATCH_SIZE * sizeof(char));

        handleErrs(
            cudaMemcpy(hostBacktrackingStringRet, deviceBacktrackStringRet, (stringLengthMax * 3) * BATCH_SIZE * sizeof(char), cudaMemcpyDeviceToHost),
            "FAILED TO COPY BACKTRACKING STRING FROM DEVICE --> HOST\n"
        );

        memalloc_time += get_time() - start_memalloc;

        uint64_t start_print = get_time();
        for (int i = sequenceIdx; i < sequenceIdx+BATCH_SIZE; i+=2) {

            uint32_t packedScore = hostSimilarityScores[i-sequenceIdx];
            s16x2 unpackedScore = unpack_s16x2(packedScore);
            int16_t scoreA = unpackedScore.A;
            int16_t scoreB = unpackedScore.B;
        
            // Backtrack matrices
            printf("%d | %d\n", i, scoreA);
            int spacingA = hostStringSpacing[i-sequenceIdx];

            printf("%s\n", hostBacktrackingStringRet + spacingA);
            printf("%s\n", hostBacktrackingStringRet + stringLengthMax + spacingA);
            printf("%s\n", hostBacktrackingStringRet + stringLengthMax + stringLengthMax + spacingA);

            printf("%d | %d\n", i+1, scoreB);
            int spacingB = hostStringSpacing[i+1-sequenceIdx];
            printf("%s\n", hostBacktrackingStringRet + spacingB);
            printf("%s\n", hostBacktrackingStringRet + stringLengthMax + spacingB);
            printf("%s\n", hostBacktrackingStringRet + stringLengthMax + stringLengthMax + spacingB);
        }
        printing_time += get_time() - start_print;

        free(hostBacktrackingStringRet);
        cudaFree(deviceBacktrackStringRet);
        cudaFree(deviceMatricesAll);
    }

    cudaFree(deviceSequences);
    cudaFree(deviceAllSequenceInfo);

    free(hostBacktrackingIndices);
    free(hostSimilarityScores);
    free(hostStringSpacing);

    cudaFree(deviceBacktrackingIndices);
    cudaFree(deviceSimilarityScores);
    cudaFree(deviceStringSpacing);


    uint64_t elapsed_time = get_elapsed_time();
    printf("Max reference length: %d\n", fileInfo.maxReferenceLength);
    printf("Max query length: %d\n", fileInfo.maxQueryLength);
    printf("Min reference length: %d\n", fileInfo.minReferenceLength);
    printf("Min query length: %d\n", fileInfo.minQueryLength);
    printf("Avg reference length: %f\n", fileInfo.avgReferenceLength);
    printf("Avg query length: %f\n", fileInfo.avgQueryLength);

    printf("Number of cells: %d\n", fileInfo.numCells);
    double kernel_time_double = ((double)kernel_time) / 1000000;
    printf("Kernel time (sec): %f\n", kernel_time_double);
    double GCUPS = (fileInfo.numCells / kernel_time_double) / 1000000000;
    printf("GCUPS: %f\n", GCUPS);

    printf("Elapsed time (usec): %lld\n", elapsed_time);
    printf("Elapsed kernel time (usec): %lld\n", kernel_time);
    printf("Elapsed backtracking time (usec): %lld\n", backtracking_time);
    printf("Elapsed memallocing time (usec): %lld\n", memalloc_time);
    printf("Elapsed printing time (usec): %lld\n", printing_time);
    printf("Elapsed time sum (usec): %lld\n",kernel_time + backtracking_time + memalloc_time + printing_time);

    // Cleanup
    printf("Cleaning up\n");
    cleanupParsedFile(allSequenceInfo, sequences);
}