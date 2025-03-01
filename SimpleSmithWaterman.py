import numpy as np
from SequenceAligner import SequenceAligner

# Performing Simple Smith Waterman
# Local Sequence Alignment

# Differences from Needleman-Wunsch --> Smith-Waterman
#   i. Initialization is all 0s instead of subtracting by gap score
#  ii. Elements CANNOT be negative (ReLU activation)
# iii. Similarity score is determine by best score within the entire matrix, not the end score
#  iv. Backtracking starts at the highest similarity score

class SmithWatermanAligner(SequenceAligner):
  
  def __init__(self, reference, query, MATCH, MISMATCH, GAP):
    super().__init__(reference, query)
    
    self.MATCH = MATCH         # Diagonal Matches
    self.MISMATCH = MISMATCH   # Diagonal Mismatches
    self.GAP = GAP             # Vertical/Horizontal Insertion/Deletions
    
  
  def execute(self):
    # Initialize Memo correctly
    # Memo size (len(self.query) + 1) * (len(self.reference) + 1)
    self.initializeMemoMatrix()

    print("\nInitialized Memo Matrix:") 
    self.printMemoMatrix()

    self.performRecursiveAnalysis()

    print("\nFinal Memo Matrix:")
    self.printMemoMatrix()
    
    # Performing backtracking to get longest subsequence
    trackerX, trackerY, trackerStar, similarityScore = self.backtrack()
    
    print()
    print("TrackerX:", trackerX)
    print("TrackerS:", trackerStar)
    print("TrackerY:", trackerY)
    print("Similarity Score:", similarityScore)


  # Memo size (len(self.query) + 1) * (len(self.reference) + 1)
  def initializeMemoMatrix(self) -> None:
    
    # NEW: All values are initialized to 0
    # Focussing on local alignment os it doesn't matter where you start
    self.Memo = np.zeros((len(self.query) + 1, len(self.reference) + 1))
    

  # Fill the 2D Memo matrix
  def performRecursiveAnalysis(self) -> None:
    
    print("\nPerforming LCS Ananlysis")
    print("Reference:", self.reference)
    print("Query:", self.query)
    print(f"(Match, Mismatch, Gap) = ({self.MATCH}, {self.MISMATCH}, {self.GAP})")
    
    for row_idx in range(1, len(self.query) + 1):
      for col_idx in range(1, len(self.reference) + 1):
        
        # Calculating max of either GAP score
        self.Memo[row_idx][col_idx] = max(
                        self.Memo[row_idx - 1][col_idx] + self.GAP, 
                        self.Memo[row_idx][col_idx - 1] + self.GAP)
        
        # If match, calculate which one is larger
        if (self.query[row_idx - 1] == self.reference[col_idx - 1]):
          match_score = self.Memo[row_idx - 1][col_idx - 1] + self.MATCH
          self.Memo[row_idx][col_idx] = max(self.Memo[row_idx][col_idx], match_score)
        
        # If mismatch, calculate which one is larger
        else:
          mismatch_score = self.Memo[row_idx - 1][col_idx - 1] + self.MISMATCH
          self.Memo[row_idx][col_idx] = max(self.Memo[row_idx][col_idx], mismatch_score)
        
        # NEW: Each element has a minimum value of 0
        #      No element in Smith-Waterman can be negative
        self.Memo[row_idx][col_idx] = max(0, self.Memo[row_idx][col_idx])
          
    return


  # Performing backtracking to get longest subsequence
  def backtrack(self):
    
    # Finding location of largest similarity on 2D array
    a = np.unravel_index(self.Memo.argmax(), self.Memo.shape)
    
    bestSimilarityScore = self.Memo[a]
    
    print()
    
    # Performing backtracking prcoess
    currentY = a[0]
    currentX = a[1]
    
    trackerX = ""
    trackerY = ""
    trackerStar = ""

    while ((currentX != 0) and (currentY != 0) and self.Memo[currentY][currentX] > 0):
      
      # print(currentY, currentX)
      
      if (self.reference[currentX-1] == self.query[currentY-1]):
        # print("MATCH")
        trackerX = self.reference[currentX-1] + trackerX
        trackerStar = "*" + trackerStar
        trackerY = self.reference[currentX-1] + trackerY
        
        currentX-=1
        currentY-=1
        
      elif (self.Memo[currentY-1][currentX-1] >= max(self.Memo[currentY][currentX-1], self.Memo[currentY-1][currentX])):
        # print("MISMATCH")
        trackerX = self.reference[currentX-1] + trackerX
        trackerStar = "|" + trackerStar
        trackerY = self.query[currentX-1] + trackerY
        
        currentX-=1
        currentY-=1
        
      elif (self.Memo[currentY][currentX-1] > self.Memo[currentY-1][currentX]):
        # print("INSERTION?")
        trackerX = "_" + trackerX
        trackerStar = " " + trackerStar
        trackerY = self.query[currentY-1] + trackerY
        
        currentY-=1
        
      else:
        print("DELETION?")
        trackerX = self.reference[currentX-1] + trackerX
        trackerStar = " " + trackerStar
        trackerY = "_" + trackerY
        
        currentX-=1
    
    return trackerX, trackerY, trackerStar, bestSimilarityScore