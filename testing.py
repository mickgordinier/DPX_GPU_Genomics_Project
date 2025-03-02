from SimpleSmithWaterman import SmithWatermanAligner
from SimpleNeedlemanWunsch import NeedlemanWunschAligner
from LongestCommonSubsequence import LongestCommonSubsequence

# sw = SmithWatermanAligner("GAATTCAGT", "GGATCGA", 5, -2, -3)
# sw.execute()

lcs = LongestCommonSubsequence("GAATTCAGT", "GGATCGA")
lcs.execute()

# nw = NeedlemanWunschAligner("GGATCGA", "GAATTCAGT", 5, -2, -3)
# nw.execute()