import numpy as np
from collections import deque
from SequenceAligner import SequenceAligner

# Find the longest subseqences between 2 strings (self.reference and self.query)
# Explination on LCS: https://www.geeksforgeeks.org/longest-common-subsequence-dp-4/
# Basically Needleman-Wunsch, but with no penalties for gaps, insertions, and deletions

class LongestCommonSubsequence(SequenceAligner):
  
  def __init__(self, reference, query):
    super().__init__(reference, query)
    
  
  def execute(self):
    
    print("\nExecuting Default Longest Common Subsequence")
    print("Reference:", self.reference)
    print("Query:", self.query)
    
    # Initialize Memo correctly
    # Memo size (len(self.query) + 1) * (len(self.reference) + 1)
    self.initializeMemoMatrix()

    self.printMemoMatrix("Initialized Memo Matrix:")

    self.performRecursiveAnalysis()

    self.printMemoMatrix("Final Memo Matrix:")
    
    # Performing backtracking and printing out all the subsequences
    self.backtrackPrintAllPaths()
    
    
  # self.Memo size (len(self.query) + 1) * (len(self.reference) + 1)
  def initializeMemoMatrix(self) -> None:
    self.Memo = np.zeros((len(self.query) + 1, len(self.reference) + 1))
  
  
  # Fill the 2D self.Memo matrix
  def performRecursiveAnalysis(self) -> None:
    
    for row_idx in range(1, len(self.query) + 1):
      for col_idx in range(1, len(self.reference) + 1):

        if self.query[row_idx - 1] == self.reference[col_idx - 1]:
          self.Memo[row_idx][col_idx] = self.Memo[row_idx - 1][col_idx - 1] + 1
        else:
          self.Memo[row_idx][col_idx] = max(self.Memo[row_idx - 1][col_idx], self.Memo[row_idx][col_idx - 1])
        
    return
  
  # Performing backtracking to get longest subsequence
  def backtrackPrintAllPaths(self) -> None:
    
    # Appending the bottom right location to the queue
    # That will always be the starting point for needleman-wunsch
    # The index is pointing to the last element of the matrix, not of the strings
    currentReferenceIdx = len(self.reference)
    currentQueryIdx = len(self.query)
    trackerLCS = ""
    
    pathsQueue = deque()
    
    pathsQueue.append([
      currentReferenceIdx, 
      currentQueryIdx, 
      trackerLCS]
    )
    
    # While there are still paths to traverse
    while (pathsQueue):
      
      path = pathsQueue.popleft()
      
      currentReferenceIdx = path[0]
      currentQueryIdx = path[1]
      trackerLCS = path[2]
      
      # If either of the indices reaches 0, then there is no more matches to be found
      if ((currentReferenceIdx != 0) and (currentQueryIdx != 0)):
        
        if (self.reference[currentReferenceIdx-1] == self.query[currentQueryIdx-1]):
          
          pathsQueue.append([
            currentReferenceIdx-1, 
            currentQueryIdx-1, 
            self.reference[currentReferenceIdx-1] + trackerLCS]
          )
          continue
        
        # DELIBERATELY 2 IF STATEMENTS
        # Both are true if they are equal --> Need to check both paths (If no match)
        if (self.Memo[currentQueryIdx][currentReferenceIdx-1] >= self.Memo[currentQueryIdx-1][currentReferenceIdx]):
          
          pathsQueue.append([
            currentReferenceIdx-1, 
            currentQueryIdx, 
            trackerLCS]
          )
          
        if (self.Memo[currentQueryIdx][currentReferenceIdx-1] <= self.Memo[currentQueryIdx-1][currentReferenceIdx]):
          
          pathsQueue.append([
            currentReferenceIdx, 
            currentQueryIdx-1, 
            trackerLCS]
          )
      
      else:
        # Given Path Trace has reached  the border and is finished
        print()
        print(trackerLCS)
    print()