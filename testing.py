from SimpleSmithWaterman import SmithWatermanAligner
from SimpleNeedlemanWunsch import NeedlemanWunschAligner
from LongestCommonSubsequence import LongestCommonSubsequence

# sw = SmithWatermanAligner("GAATTCAGT", "GGATCGA", 5, -2, -3)
# sw.execute()

# lcs = LongestCommonSubsequence("GAATTCAGT", "GGATCGA")
# lcs.execute()

nw = NeedlemanWunschAligner("GGATCGA", "GAATTCAGT", 5, -2, -3)
nw.execute()

# Initialize self.Memo correctly
# self.Memo size (len(self.query) + 1) * (len(self.reference) + 1)
# nw.initializeMemoMatrix()

# print("\nInitialized self.Memo Matrix:") 
# nw.printMemoMatrix()

# nw.performRecursiveAnalysis()

# print("\nFinal self.Memo Matrix:")
# nw.printMemoMatrix()

# # Performing backtracking to get longest subsequence
# nw.backtrack()

# print()
# print("TrackerX:", trackerX)
# print("TrackerY:", trackerY)
# print("Similarity Score:", nw.Memo[len(nw.query)][len(nw.reference)])
