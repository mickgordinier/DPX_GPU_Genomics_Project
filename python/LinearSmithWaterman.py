import numpy as np
from collections import deque
from SequenceAligner import SequenceAligner

# Performing Simple Smith Waterman
# Local Sequence Alignment

CORNER = 0
LEFT_GAP = 1
UPPER_GAP = 2

# Differences from Needleman-Wunsch --> Smith-Waterman
#   i. Initialization is all 0s instead of subtracting by gap score
#  ii. Elements CANNOT be negative (ReLU activation)
# iii. Similarity score is determine by best score within the entire matrix, not the end score
#  iv. Backtracking starts at the highest similarity score

class LinearSmithWatermanAligner(SequenceAligner):
  
  def __init__(self, reference, query, MATCH, MISMATCH, GAP):
    super().__init__(reference, query)
    
    self.MATCH = MATCH         # Diagonal Matches
    self.MISMATCH = MISMATCH   # Diagonal Mismatches
    self.GAP = GAP             # Vertical/Horizontal Insertion/Deletions
    
  
  def execute(self):
    
    print("\nExecuting Linear Gap Penalty Smith-Waterman")
    print("Reference:", self.reference)
    print("Query:", self.query)
    print(f"(Match, Mismatch, Gap) = ({self.MATCH}, {self.MISMATCH}, {self.GAP})")
    
    # Initialize Memo correctly
    # Memo size (len(self.query) + 1) * (len(self.reference) + 1)
    self.initializeMemoMatrix()

    self.printMemoMatrix("Initialized Memo Matrix:")

    self.performRecursiveAnalysis()

    self.printMemoMatrix("Final Memo Matrix:")
    
    # Performing backtracking and printing out all the subsequences
    self.backtrackPrintAllPaths()


  # Memo size (len(self.query) + 1) * (len(self.reference) + 1)
  def initializeMemoMatrix(self) -> None:
    
    # NEW: All values are initialized to 0
    # Focussing on local alignment os it doesn't matter where you start
    self.Memo = np.zeros((len(self.query) + 1, len(self.reference) + 1))
    

  # Fill the 2D Memo matrix
  def performRecursiveAnalysis(self) -> None:
    
    # Building a backtracking matrix to keep track of where it came from
    # Each element has 3 places it could backtrack to [corner, left, up]
    # Keeping track for each element whether or not it came from that location
    # It is possible to come from multiple locations (branches of backtracking)
    self.backtrackMatrix = np.zeros((len(self.query), len(self.reference), 3), dtype=bool)
    
    for row_idx in range(1, len(self.query) + 1):
      for col_idx in range(1, len(self.reference) + 1):
        
        upper_gap = self.Memo[row_idx - 1][col_idx] + self.GAP
        left_gap = self.Memo[row_idx][col_idx - 1] + self.GAP
        
        # If match, add the match score from the corner
        if (self.query[row_idx - 1] == self.reference[col_idx - 1]):
          corner_score = self.Memo[row_idx - 1][col_idx - 1] + self.MATCH
        # Otherwise, add the penalty score from the corner
        else:
          corner_score = self.Memo[row_idx - 1][col_idx - 1] + self.MISMATCH
          
        temp_greatest = max(upper_gap, left_gap, corner_score)
        
        # NEW: Each element has a minimum value of 0
        #      No element in Smith-Waterman can be negative
        self.Memo[row_idx][col_idx] = max(0, temp_greatest)
        
        # If the greatest amongst the 3 is still negative, no backtracking
        if (temp_greatest < 0):
          self.backtrackMatrix[row_idx-1][col_idx-1][CORNER] = 0
          self.backtrackMatrix[row_idx-1][col_idx-1][LEFT_GAP] = 0
          self.backtrackMatrix[row_idx-1][col_idx-1][UPPER_GAP] = 0
        
        # Determining where max value came from
        else:
          self.backtrackMatrix[row_idx-1][col_idx-1][CORNER] = (corner_score == self.Memo[row_idx][col_idx])
          self.backtrackMatrix[row_idx-1][col_idx-1][LEFT_GAP] = (left_gap == self.Memo[row_idx][col_idx])
          self.backtrackMatrix[row_idx-1][col_idx-1][UPPER_GAP] = (upper_gap == self.Memo[row_idx][col_idx]) 
        # end if
        
      # end for col_idx
    # end for row_idx
    
    return


  # Performing backtracking to get longest subsequence
  def backtrackPrintAllPaths(self):
    
    # Finding location of largest similarity on 2D array
    
    bestSimilarityScore = np.max(self.Memo)
    rowIndices, colIndices = np.where(self.Memo == bestSimilarityScore)
    
    # Appending the bottom right location to the queue
    # That will always be the starting point for needleman-wunsch
    # The index is pointing to the last element of the matrix, not of the strings
    pathsQueue = deque()
    
    for max_score_idx in range(len(rowIndices)):
      
      currentReferenceIdx = colIndices[max_score_idx]
      currentQueryIdx = rowIndices[max_score_idx]
      trackerReference = ""
      trackerConnection = ""
      trackerQuery = ""
      
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
      if ((currentQueryIdx != 0) and
          (currentReferenceIdx != 0) and 
          (self.Memo[currentQueryIdx][currentReferenceIdx] > 0)):
        
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
  