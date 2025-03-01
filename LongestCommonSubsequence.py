import numpy as np
from SequenceAligner import SequenceAligner

# Find the longest subseqences between 2 strings (self.reference and self.query)
# Explination on LCS: https://www.geeksforgeeks.org/longest-common-subsequence-dp-4/
# Basically Needleman-Wunsch, but with no penalties for gaps, insertions, and deletions

class LongestCommonSubsequence(SequenceAligner):
  
  def __init__(self, reference, query):
    super().__init__(reference, query)
    
    
  # Performs Longest Common Substring on 2 Strings
  def execute(self) -> None:

    # Initialize self.Memo correctly
    # self.Memo size (len(self.query) + 1) * (len(self.reference) + 1)
    self.initializeMemoMatrix()

    print("\nInitialized self.Memo Matrix:") 
    self.printMemoMatrix()

    self.performRecursiveAnalysis()

    print("\nFinal self.Memo Matrix:")
    self.printMemoMatrix()
    
    # Performing backtracking to get longest subsequence
    lcs = self.backtrack()
    
    print()
    print("Lonest Subsequence:", lcs)
    
    
  # self.Memo size (len(self.query) + 1) * (len(self.reference) + 1)
  def initializeMemoMatrix(self) -> None:
    self.Memo = np.zeros((len(self.query) + 1, len(self.reference) + 1))
  
  
  # Fill the 2D self.Memo matrix
  def performRecursiveAnalysis(self) -> None:
    
    print("\nPerforming LCS Ananlysis")
    print("Reference:", self.reference)
    print("Query:", self.query)
    
    for row_idx in range(1, len(self.query) + 1):
      for col_idx in range(1, len(self.reference) + 1):

        if self.query[row_idx - 1] == self.reference[col_idx - 1]:
          self.Memo[row_idx][col_idx] = self.Memo[row_idx - 1][col_idx - 1] + 1
        else:
          self.Memo[row_idx][col_idx] = max(self.Memo[row_idx - 1][col_idx], self.Memo[row_idx][col_idx - 1])
        
    return
  
  # Performing backtracking to get longest subsequence
  def backtrack(self) -> None:
    
    currentX = len(self.reference)
    currentY = len(self.query)

    tracker = ""

    while ((currentX != 0) and (currentY != 0)):
      if (self.reference[currentX-1] == self.query[currentY-1]):
        tracker = self.reference[currentX-1] + tracker
        currentX-=1
        currentY-=1
      elif (self.Memo[currentY][currentX-1] > self.Memo[currentY-1][currentX]):
        currentX-=1
      else:
        currentY-=1
    
    return tracker