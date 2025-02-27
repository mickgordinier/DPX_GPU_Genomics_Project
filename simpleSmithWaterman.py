import numpy as np

# Performing Simple Smith Waterman
# Local Sequence Alignment

# Differences from Needleman-Wunsch --> Smith-Waterman
#   i. Initialization is all 0s instead of subtracting by gap score
#  ii. Elements CANNOT be negative (ReLU activation)
# iii. Similarity score is determine by best score within the entire matrix, not the end score
#  iv. Backtracking starts at the highest similarity score

X = "GAATTCAGT"
Y = "GGATCGA"

# Assigning scores for each scenarios
# Assuming scores to be constant throughout (for simplicity)
MATCH = 5       # Diagonal Matches
MISMATCH = -2   # Diagonal Mismatches
GAP = -3        # Vertical/Horizontal Insertion/Deletions


# Prints the matrix in a nice format
def printLccMatrix(Memo: np.ndarray[np.any], X: str, Y: str) -> str:
  
  lineToPrint = "      "
  for char in X:
      lineToPrint += char + "  "
  print(lineToPrint)

  print(" ", Memo[0])
  for abc in range(1, len(Y) + 1):
    print(Y[abc - 1], Memo[abc])
  
  return


# Memo size (len(Y) + 1) * (len(X) + 1)
def initializeMemoMatrix() -> np.ndarray[np.any]:
  
  Memo = np.zeros((len(Y) + 1, len(X) + 1))
  
  # NEW: All values are initialized to 0
  # Focussing on local alignment os it doesn't matter where you start
    
  return Memo


# Fill the 2D Memo matrix
def populateLcmMatrix(Memo: np.ndarray[np.any], X: str, Y: str) -> None:
  
  for row_idx in range(1, len(Y) + 1):
    for col_idx in range(1, len(X) + 1):
      
      # Calculating max of either GAP score
      Memo[row_idx][col_idx] = max(Memo[row_idx - 1][col_idx] + GAP, Memo[row_idx][col_idx - 1] + GAP)
      
      # If match, calculate which one is larger
      if (Y[row_idx - 1] == X[col_idx - 1]):
        match_score = Memo[row_idx - 1][col_idx - 1] + MATCH
        Memo[row_idx][col_idx] = max(Memo[row_idx][col_idx], match_score)
      
      # If mismatch, calculate which one is larger
      else:
        mismatch_score = Memo[row_idx - 1][col_idx - 1] + MISMATCH
        Memo[row_idx][col_idx] = max(Memo[row_idx][col_idx], mismatch_score)
      
      # NEW: Each element has a minimum value of 0
      #      No element in Smith-Waterman can be negative
      Memo[row_idx][col_idx] = max(0, Memo[row_idx][col_idx])
        
  return


# Performing backtracking to get longest subsequence
def backtrack(Memo: np.ndarray[np.any], X: str, Y: str):
  
  # Finding location of largest similarity on 2D array
  a = np.unravel_index(Memo.argmax(), Memo.shape)
  
  bestSimilarityScore = Memo[a]
  
  print()
  
  # Performing backtracking prcoess
  currentY = a[0]
  currentX = a[1]
  
  trackerX = ""
  trackerY = ""
  trackerStar = ""

  while ((currentX != 0) and (currentY != 0) and Memo[currentY][currentX] > 0):
    
    print(currentY, currentX)
    
    if (X[currentX-1] == Y[currentY-1]):
      print("MATCH")
      trackerX = X[currentX-1] + trackerX
      trackerStar = "*" + trackerStar
      trackerY = X[currentX-1] + trackerY
      
      currentX-=1
      currentY-=1
      
    elif (Memo[currentY-1][currentX-1] >= max(Memo[currentY][currentX-1], Memo[currentY-1][currentX])):
      print("MISMATCH")
      trackerX = X[currentX-1] + trackerX
      trackerStar = "|" + trackerStar
      trackerY = Y[currentX-1] + trackerY
      
      currentX-=1
      currentY-=1
      
    elif (Memo[currentY][currentX-1] > Memo[currentY-1][currentX]):
      print("INSERTION?")
      trackerX = "_" + trackerX
      trackerStar = " " + trackerStar
      trackerY = Y[currentY-1] + trackerY
      
      currentY-=1
      
    else:
      print("DELETION?")
      trackerX = X[currentX-1] + trackerX
      trackerStar = " " + trackerStar
      trackerY = "_" + trackerY
      
      currentX-=1
  
  return trackerX, trackerY, trackerStar, bestSimilarityScore


# Performs Longest Common Substring on 2 Strings
def performSmithWaterman(X: str, Y:str) -> None:

  # Initialize Memo correctly
  # Memo size (len(Y) + 1) * (len(X) + 1)
  Memo = initializeMemoMatrix()

  # Print header for X
  print("\nInitialized Memo Matrix:") 
  printLccMatrix(Memo, X, Y)

  # Fill the Memo matrix
  populateLcmMatrix(Memo, X, Y)

  # Final print of the memoized matrix after it's fully populated
  print("\nFinal Memo Matrix:")
  printLccMatrix(Memo, X, Y)
  
  # Performing backtracking to get longest subsequence
  trackerX, trackerY, trackerStar, similarityScore = backtrack(Memo, X, Y)
  
  print()
  print("TrackerX:", trackerX)
  print("TrackerS:", trackerStar)
  print("TrackerY:", trackerY)
  print("Similarity Score:", similarityScore)


performSmithWaterman(X, Y)