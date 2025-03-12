import numpy as np
from collections import deque
from SequenceAligner import SequenceAligner

# Performing Simple Needleman Wunsch
# Global Sequence Alignment

RED = '\033[31m'
GREEN = '\033[32m'
RESET = '\033[0m'

CORNER = 0
LEFT_GAP = 1
UPPER_GAP = 2

# Differences from Longest Common Subsequence --> Needleman-Wunsch
#   i. Now penalties are enforced (GAP and MISMATCH)
#  ii. Inclusion of different MATCH, MISMATCH, and GAP score

class AffineNeedlemanWunschAligner(SequenceAligner):
  
  def __init__(self, reference, query, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND):
    super().__init__(reference, query)
    
    self.MATCH = MATCH         # Diagonal Matches
    self.MISMATCH = MISMATCH   # Diagonal Mismatches
    self.GAP_OPEN = GAP_OPEN      
    self.GAP_EXTEND = GAP_EXTEND  
  
  
  # Prints the matrix in a nice format
  def printMemoMatrix(self, title:str) -> None:
    
    print("\n============================================")
    print(title, RED + "REFERENCE " + GREEN + "QUERY" + RESET)
    print("--------------------------------------------\n")
    
    lineToPrint = "            "
    for char in self.reference:
        lineToPrint += char + "     "
    
    
    # ----- PRINTING DISTANCE MATRIX ---- #
    
    print("MATRIX DISTANCE")
    print(RED + lineToPrint + RESET)
      
    # Print numpy matrix using the numpy matrix print and save it to a string
    # Pad each digit to 3 decimal places for nice formatting
    matrix_str = np.array2string(self.Distance, 
                             formatter={'all': lambda x: f'{int(x):4d}' if np.isfinite(x) else str(x)}, 
                             separator=', ', 
                             max_line_width=80, 
                             precision=0)
    
    # Remove the beginning and end backets
    matrix_str = matrix_str[1:-1]
    
    # Remove the spacing after each newline representing a new row
    matrix_rows = [row.strip() for row in matrix_str.split("\n")]
    
    # Add each character of query to each row
    matrix_rows[0] = "  " + matrix_rows[0]
    for row_idx in range(1, len(self.query)+1):
      matrix_rows[row_idx] = f"{GREEN}{self.query[row_idx-1]}{RESET} {matrix_rows[row_idx]}"
    
    matrix_str = "\n".join(matrix_rows)
    
    print(matrix_str)
    print()
    
    # ----- PRINTING Best_Reference_Gaps_End MATRIX ---- #
    
    print("MATRIX Best_Reference_Gaps_End")
    print(RED + lineToPrint + RESET)
      
    # Print numpy matrix using the numpy matrix print and save it to a string
    # Pad each digit to 3 decimal places for nice formatting
    matrix_str = np.array2string(self.Best_Reference_Gaps_End, 
                             formatter={'all': lambda x: f'{int(x):4d}' if np.isfinite(x) else str(x)}, 
                             separator=', ', 
                             max_line_width=80, 
                             precision=0)
    
    # Remove the beginning and end backets
    matrix_str = matrix_str[1:-1]
    
    # Remove the spacing after each newline representing a new row
    matrix_rows = [row.strip() for row in matrix_str.split("\n")]
    
    # Add each character of query to each row
    matrix_rows[0] = "  " + matrix_rows[0]
    for row_idx in range(1, len(self.query)+1):
      matrix_rows[row_idx] = f"{GREEN}{self.query[row_idx-1]}{RESET} {matrix_rows[row_idx]}"
    
    matrix_str = "\n".join(matrix_rows)
    
    print(matrix_str)
    print()
    
    # ----- PRINTING Best_Query_Gaps_End MATRIX ---- #
    
    print("MATRIX Best_Query_Gaps_End")
    print(RED + lineToPrint + RESET)
      
    # Print numpy matrix using the numpy matrix print and save it to a string
    # Pad each digit to 3 decimal places for nice formatting
    matrix_str = np.array2string(self.Best_Query_Gaps_End, 
                             formatter={'all': lambda x: f'{int(x):4d}' if np.isfinite(x) else str(x)}, 
                             separator=', ', 
                             max_line_width=80, 
                             precision=0)
    
    # Remove the beginning and end backets
    matrix_str = matrix_str[1:-1]
    
    # Remove the spacing after each newline representing a new row
    matrix_rows = [row.strip() for row in matrix_str.split("\n")]
    
    # Add each character of query to each row
    matrix_rows[0] = "  " + matrix_rows[0]
    for row_idx in range(1, len(self.query)+1):
      matrix_rows[row_idx] = f"{GREEN}{self.query[row_idx-1]}{RESET} {matrix_rows[row_idx]}"
    
    matrix_str = "\n".join(matrix_rows)
    
    print(matrix_str)
    
    print("\n============================================\n")
        
    
  
  def execute(self):
    
    print("\nExecuting Default Needleman-Wunsch")
    print("Reference:", self.reference)
    print("Query:", self.query)
    print(f"(Match, Mismatch, Gap Open, Gap Extend) = ({self.MATCH}, {self.MISMATCH}, {self.GAP_OPEN}, {self.GAP_EXTEND})")
    
    # Initialize self.Memo correctly
    self.initializeMemoMatrix()

    self.printMemoMatrix("Initialized Memo Matrix:")

    self.performRecursiveAnalysis()

    self.printMemoMatrix("Final Memo Matrix:")
    
    # # Performing backtracking and printing out all the subsequences
    # self.backtrackPrintAllPaths()


  def initializeMemoMatrix(self) -> None:
    
    # The main (middle) distance matrix is of size (N+1)(M+1)
    self.Distance = np.zeros((len(self.query) + 1, len(self.reference) + 1))
    
    # Initialize the 0th row and 0th column to be adding on a penalty extension
    # Decided to add both extension and open penalty to first open
    for reference_idx in range(1, len(self.reference)+1):
      self.Distance[0][reference_idx] = ((reference_idx) * self.GAP_EXTEND) + self.GAP_OPEN
      
    for query_idx in range(1, len(self.query)+1):
      self.Distance[query_idx][0] = ((query_idx) * self.GAP_EXTEND) + self.GAP_OPEN
    
    # The best score under the additional constraint that the alignment ends in a gap within REFERENCE sequence
    self.Best_Reference_Gaps_End = np.zeros((len(self.query) + 1, len(self.reference) + 1))
    for reference_idx in range(1, len(self.reference)+1):
      self.Best_Reference_Gaps_End[0][reference_idx] = -np.inf
    
    # The best score under the additional constraint that the alignment ends in a gap within QUERY sequence
    self.Best_Query_Gaps_End = np.zeros((len(self.query) + 1, len(self.reference) + 1))
    for query_idx in range(1, len(self.query)+1):
      self.Best_Query_Gaps_End[query_idx][0] = -np.inf


  # Fill the 2D self.Memo matrix
  def performRecursiveAnalysis(self) -> None:
    
    # Building a backtracking matrix to keep track of where it came from
    # Each element has 3 places it could backtrack to [corner, left, up]
    # Keeping track for each element whether or not it came from that location
    # It is possible to come from multiple locations (branches of backtracking)
    # self.backtrackMatrix = np.zeros((len(self.query), len(self.reference), 3), dtype=bool)
  
    for query_idx in range(1, len(self.query) + 1):
      for reference_idx in range(1, len(self.reference) + 1):
        
        # The best score under the additional constraint that the alignment ends in a gap within QUERY sequence
        self.Best_Reference_Gaps_End[query_idx][reference_idx] = max(
          self.Distance[query_idx-1][reference_idx] + self.GAP_OPEN + self.GAP_EXTEND,
          self.Best_Reference_Gaps_End[query_idx-1][reference_idx] + self.GAP_EXTEND
        )
        
        # The best score under the additional constraint that the alignment ends in a gap within REFERENCE sequence
        self.Best_Query_Gaps_End[query_idx][reference_idx] = max(
          self.Distance[query_idx][reference_idx-1] + self.GAP_OPEN + self.GAP_EXTEND,
          self.Best_Query_Gaps_End[query_idx][reference_idx-1] + self.GAP_EXTEND
        )
        
        if (self.query[query_idx-1] == self.reference[reference_idx-1]):
          match_score = self.MATCH
        else:
          match_score = self.MISMATCH
        
        self.Distance[query_idx][reference_idx] = max(
          self.Distance[query_idx-1][reference_idx-1] + match_score,
          self.Best_Reference_Gaps_End[query_idx][reference_idx],
          self.Best_Query_Gaps_End[query_idx][reference_idx]
        )
        
        
        # upper_gap = self.Memo[row_idx - 1][col_idx] + self.GAP
        # left_gap = self.Memo[row_idx][col_idx - 1] + self.GAP
        
        # # If match, add the match score from the corner
        # if (self.query[row_idx - 1] == self.reference[col_idx - 1]):
        #   corner_score = self.Memo[row_idx - 1][col_idx - 1] + self.MATCH
        # # Otherwise, add the penalty score from the corner
        # else:
        #   corner_score = self.Memo[row_idx - 1][col_idx - 1] + self.MISMATCH
        
        # self.Memo[row_idx][col_idx] = max(upper_gap, left_gap, corner_score)
        
        # self.backtrackMatrix[row_idx-1][col_idx-1][CORNER] = (corner_score == self.Memo[row_idx][col_idx])
        # self.backtrackMatrix[row_idx-1][col_idx-1][LEFT_GAP] = (left_gap == self.Memo[row_idx][col_idx])
        # self.backtrackMatrix[row_idx-1][col_idx-1][UPPER_GAP] = (upper_gap == self.Memo[row_idx][col_idx])  
      # end for col_idx
    # end for row_idx
    
    return
    

  # Performing backtracking to get longest subsequence
  def backtrackPrintAllPaths(self) -> None:
    
    # Appending the bottom right location to the queue
    # That will always be the starting point for needleman-wunsch
    # The index is pointing to the last element of the matrix, not of the strings
    currentReferenceIdx = len(self.reference)
    currentQueryIdx = len(self.query)
    trackerReference = ""
    trackerConnection = ""
    trackerQuery = ""
    
    pathsQueue = deque()
    
    pathsQueue.append([
      currentReferenceIdx, 
      currentQueryIdx, 
      trackerReference, 
      trackerConnection, 
      trackerQuery]
    )
    
    # While there are still paths to traverse
    while (pathsQueue):
      
      path = pathsQueue.popleft()
      
      currentReferenceIdx = path[0]
      currentQueryIdx = path[1]
      trackerReference = path[2]
      trackerConnection = path[3]
      trackerQuery = path[4]
      
      # If either of the indices are 0, we have traversed through the sequence
      if (not((currentReferenceIdx == 0) and (currentQueryIdx == 0))):
        
        if (self.backtrackMatrix[currentQueryIdx-1][currentReferenceIdx-1][CORNER] and
            self.reference[currentReferenceIdx-1] == self.query[currentQueryIdx-1]):
          
          pathsQueue.append([
            currentReferenceIdx-1, 
            currentQueryIdx-1, 
            self.reference[currentReferenceIdx-1] + trackerReference, 
            "*" + trackerConnection, 
            self.query[currentQueryIdx-1] + trackerQuery]
          )
        #end if match
          
        if (self.backtrackMatrix[currentQueryIdx-1][currentReferenceIdx-1][CORNER] and
            self.reference[currentReferenceIdx-1] != self.query[currentQueryIdx-1]):
          
          pathsQueue.append([
            currentReferenceIdx-1, 
            currentQueryIdx-1, 
            self.reference[currentReferenceIdx-1] + trackerReference, 
            "|" + trackerConnection, 
            self.query[currentQueryIdx-1] + trackerQuery]
          )
        #end if mismatch
          
        if (self.backtrackMatrix[currentQueryIdx-1][currentReferenceIdx-1][LEFT_GAP]):
          
          pathsQueue.append([
            currentReferenceIdx-1, 
            currentQueryIdx, 
            self.reference[currentReferenceIdx-1] + trackerReference, 
            " " + trackerConnection, 
            "_" + trackerQuery]
          )
        #end if left gap
          
        if (self.backtrackMatrix[currentQueryIdx-1][currentReferenceIdx-1][UPPER_GAP]):
          
          pathsQueue.append([
            currentReferenceIdx, 
            currentQueryIdx-1, 
            "_" + trackerReference, 
            " " + trackerConnection, 
            self.query[currentQueryIdx-1] + trackerQuery]
          )
        #end if upper gap
      
      # end if trace not completed
          
      else:
        # Given Path Trace has reached  the border and is finished
        print()
        print(trackerReference)
        print(trackerConnection)
        print(trackerQuery)
        
    # end while
    print()
    
    return
  
  
anw = AffineNeedlemanWunschAligner("ABxxxCDE", "ABCDE", 7, -2, -10, -5)
anw.execute()