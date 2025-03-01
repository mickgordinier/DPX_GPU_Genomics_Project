from simpleSmithWaterman import SmithWatermanAligner
from LongestCommonSubsequence import LongestCommonSubsequence

# sw = SmithWatermanAligner("GAATTCAGT", "GGATCGA", 5, -2, -3)
# sw.execute()

lcs = LongestCommonSubsequence("GAATTCAGT", "GGATCGA")
lcs.execute()
