import numpy as np
from SequenceAligner import SequenceAligner

# Performing Simple Needleman Wunsch
# Global Sequence Alignment

# Differences from Longest Common Subsequence --> Needleman-Wunsch
#   i. Now penalties are enforced (GAP and MISMATCH)
#  ii. Inclusion of different MATCH, MISMATCH, and GAP score

class NeedlemanWunschAligner(SequenceAligner):
  
  def __init__(self, reference, query, MATCH, MISMATCH, GAP):
    super().__init__(reference, query)
    
    self.MATCH = MATCH         # Diagonal Matches
    self.MISMATCH = MISMATCH   # Diagonal Mismatches
    self.GAP = GAP             # Vertical/Horizontal Insertion/Deletions
    
  
  def execute(self):
    
    # Initialize self.Memo correctly
    # self.Memo size (len(self.query) + 1) * (len(self.reference) + 1)
    self.initializeMemoMatrix()

    print("\nInitialized self.Memo Matrix:") 
    self.printMemoMatrix()

    self.performRecursiveAnalysis()

    print("\nFinal self.Memo Matrix:")
    self.printMemoMatrix()
    
    # Performing backtracking to get longest subsequence
    trackerX, trackerY = self.backtrack()
    
    print()
    print("TrackerX:", trackerX)
    print("TrackerY:", trackerY)
    print("Similarity Score:", self.Memo[len(self.query)][len(self.reference)])


  # self.Memo size (len(self.query) + 1) * (len(self.reference) + 1)
  def initializeMemoMatrix(self) -> None:
    
    self.Memo = np.zeros((len(self.query) + 1, len(self.reference) + 1))
    
    for x_idx in range(0, len(self.reference)+1):
      self.Memo[0][x_idx] = x_idx * self.GAP
      
    for y_idx in range(0, len(self.query)+1):
      self.Memo[y_idx][0] = y_idx * self.GAP
      

  # Fill the 2D self.Memo matrix
  def performRecursiveAnalysis(self) -> None:
    
    print("\nPerforming LCS Ananlysis")
    print("Reference:", self.reference)
    print("Query:", self.query)
    print(f"(Match, Mismatch, Gap) = ({self.MATCH}, {self.MISMATCH}, {self.GAP})")
  
    for row_idx in range(1, len(self.query) + 1):
      for col_idx in range(1, len(self.reference) + 1):
        
        # Calculating max of either GAP score
        self.Memo[row_idx][col_idx] = max(self.Memo[row_idx - 1][col_idx] + self.GAP, self.Memo[row_idx][col_idx - 1] + self.GAP)
        
        # If match, calculate which one is larger
        if (self.query[row_idx - 1] == self.reference[col_idx - 1]):
          match_score = self.Memo[row_idx - 1][col_idx - 1] + self.MATCH
          self.Memo[row_idx][col_idx] = max(self.Memo[row_idx][col_idx], match_score)
        
        # If mismatch, calculate which one is larger
        else:
          mismatch_score = self.Memo[row_idx - 1][col_idx - 1] + self.MISMATCH
          self.Memo[row_idx][col_idx] = max(self.Memo[row_idx][col_idx], mismatch_score)
        

  # Performing backtracking to get longest subsequence
  def backtrack(self):
  
    currentX = len(self.reference)
    currentY = len(self.query)

    tracker = ""
    
    trackerX = ""
    trackerY = ""

    while ((currentX != 0) and (currentY != 0)):
      
      if (self.reference[currentX-1] == self.query[currentY-1]):
        trackerX = self.reference[currentX-1] + trackerX
        trackerY = self.reference[currentX-1] + trackerY
        
        currentX-=1
        currentY-=1
        
      elif (self.Memo[currentY-1][currentX-1] >= max(self.Memo[currentY][currentX-1], self.Memo[currentY-1][currentX])):
        trackerX = self.reference[currentX-1] + trackerX
        trackerY = self.reference[currentX-1] + trackerY
        
        currentX-=1
        currentY-=1
        
      elif (self.Memo[currentY][currentX-1] > self.Memo[currentY-1][currentX]):
        trackerX = "_" + trackerX
        trackerY = self.query[currentY-1] + trackerY
        currentX-=1
        
      else:
        trackerX = self.reference[currentX-1] + trackerX
        trackerY = "_" + trackerY
        currentY-=1
    
    return trackerX, trackerY