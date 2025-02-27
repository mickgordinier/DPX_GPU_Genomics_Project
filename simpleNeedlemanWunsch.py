# Performing Simple Needleman Wunsch
# Global Sequence Alignment

# Differences from Longest Common Subsequence --> Needleman-Wunsch
#   i. Now penalties are enforced (GAP and MISMATCH)
#  ii. Inclusion of different MATCH, MISMATCH, and GAP score

X = "AACG"
Y = "AATCG"

# Assigning scores for each scenarios
# Assuming scores to be constant throughout (for simplicity)
MATCH    = 5    # Diagonal Matches
MISMATCH = -1   # Diagonal Mismatches
GAP      = -6   # Vertical/Horizontal Insertion/Deletions


# Prints the matrix in a nice format
def printLccMatrix(Memo: list[list[str]], X: str, Y: str) -> str:
  
  lineToPrint = "      "
  for char in X:
      lineToPrint += char + "  "
  print(lineToPrint)

  print(" ", Memo[0])
  for abc in range(1, len(Y) + 1):
    print(Y[abc - 1], Memo[abc])
  
  return


# Memo size (len(Y) + 1) * (len(X) + 1)
def initializeMemoMatrix() -> list[list[str]]:
  Memo = [[0] * (len(X) + 1) for _ in range(len(Y) + 1)]
  
  for x_idx in range(0, len(X)+1):
    Memo[0][x_idx] = x_idx * GAP
    
  for y_idx in range(0, len(Y)+1):
    Memo[y_idx][0] = y_idx * GAP
    
  return Memo


# Fill the 2D Memo matrix
def populateLcmMatrix(Memo: list[list[str]], X: str, Y: str) -> None:
  
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
        
  return


# Performing backtracking to get longest subsequence
def backtrack(Memo: list[list[str]], X: str, Y: str):
  
  currentX = len(X)
  currentY = len(Y)

  tracker = ""
  
  trackerX = ""
  trackerY = ""

  while ((currentX != 0) and (currentY != 0)):
    
    if (X[currentX-1] == Y[currentY-1]):
      trackerX = X[currentX-1] + trackerX
      trackerY = X[currentX-1] + trackerY
      
      currentX-=1
      currentY-=1
      
    elif (Memo[currentY-1][currentX-1] >= max(Memo[currentY][currentX-1], Memo[currentY-1][currentX])):
      trackerX = X[currentX-1] + trackerX
      trackerY = X[currentX-1] + trackerY
      
      currentX-=1
      currentY-=1
      
    elif (Memo[currentY][currentX-1] > Memo[currentY-1][currentX]):
      trackerX = "_" + trackerX
      trackerY = Y[currentY-1] + trackerY
      currentX-=1
      
    else:
      trackerX = X[currentX-1] + trackerX
      trackerY = "_" + trackerY
      currentY-=1
  
  return trackerX, trackerY


# Performs Longest Common Substring on 2 Strings
def performNeedlemanWunsch(X: str, Y:str) -> None:

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
  trackerX, trackerY = backtrack(Memo, X, Y)
  
  print()
  print("TrackerX:", trackerX)
  print("TrackerY:", trackerY)
  print("Similarity Score:", Memo[len(Y)][len(X)])


performNeedlemanWunsch(X, Y)