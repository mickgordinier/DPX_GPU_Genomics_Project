# Find the longest subseqences between 2 strings (X and Y)
# Explination on LCS: https://www.geeksforgeeks.org/longest-common-subsequence-dp-4/
# Basically Needleman-Wunsch, but with no penalties for gaps, insertions, and deletions

X = "AACG"
Y = "AATCG"

# X = "ZXCVBNMASDFGHJKLQWERTYUIOPLKJHGFDSAZXCVBNMASDFGHJKLQWERTYUIOP"
# Y = "ASDFGHJKLZXCVBNMASDFGHJKLQWERTYUIOZXCVBNMASDFGHJKLQWERTYUIOP"


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
def initializeMemoMatrix() -> None:
  return [[0] * (len(X) + 1) for _ in range(len(Y) + 1)]


# Fill the 2D Memo matrix
def populateLcmMatrix(Memo: list[list[str]], X: str, Y: str) -> None:
  
  for row_idx in range(1, len(Y) + 1):
    for col_idx in range(1, len(X) + 1):

      if Y[row_idx - 1] == X[col_idx - 1]:
        Memo[row_idx][col_idx] = Memo[row_idx - 1][col_idx - 1] + 1
      else:
        Memo[row_idx][col_idx] = max(Memo[row_idx - 1][col_idx], Memo[row_idx][col_idx - 1])
        
  return


# Performing backtracking to get longest subsequence
def backtrack(Memo: list[list[str]], X: str, Y: str) -> str:
  
  currentX = len(X)
  currentY = len(Y)

  tracker = ""

  while ((currentX != 0) and (currentY != 0)):
    if (X[currentX-1] == Y[currentY-1]):
      tracker = X[currentX-1] + tracker
      currentX-=1
      currentY-=1
    elif (Memo[currentY][currentX-1] > Memo[currentY-1][currentX]):
      currentX-=1
    else:
      currentY-=1
  
  return tracker


# Performs Longest Common Substring on 2 Strings
def performLCM(X: str, Y:str) -> None:

  # Initialize Memo correctly
  # Memo size (len(Y) + 1) * (len(X) + 1)
  Memo = initializeMemoMatrix()

  # Print header for X
  printLccMatrix(Memo, X, Y)

  # Fill the Memo matrix
  populateLcmMatrix(Memo, X, Y)

  # Final print of the memoized matrix after it's fully populated
  print("\nFinal Memo Matrix:")
  printLccMatrix(Memo, X, Y)
  
  # Performing backtracking to get longest subsequence
  lcm = backtrack(Memo, X, Y)
  
  print()
  print("Lonest Subsequence:", lcm)


performLCM(X, Y)