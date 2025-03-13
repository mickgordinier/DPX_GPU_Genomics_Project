from LinearSmithWaterman import LinearSmithWatermanAligner
from LinearNeedlemanWunsch import LinearNeedlemanWunschAligner
from LongestCommonSubsequence import LongestCommonSubsequence
from AffineNeedlemanWunsch import AffineNeedlemanWunschAligner

# This test is to showcase the close differences/similarities
# between the sequence alignment algorithms: 
#   i. Longest Common Subsequence
#  ii. Needleman-Wunsch (Global Sequence Aligner)
# iii. Smith-Waterman (Local Sequence Aligner)

#   i. Longest Common Subsequence (LCS):
# This algorithm can be thought of as a specific case of Needleman-Wunsch.
# This algorithm is a global sequence aligner, with no penalty scores.
# There is no "scoring", just wanting to MAXIMIZE LETTER COUNT (Matching)
# Ex. LCS("AxBxCxDx", "yAyyByCyDy") = "ABCD"
lcs = LongestCommonSubsequence("AxBxCxDx", "yAyyByCyDy")
lcs.execute()

#   ii. Needleman-Wunsch (Global Sequence Aligner)
# This algorithm focusses on globally aligning 2 genomic sequences ("Reference" and "Query" genomes).
# With NW, usually both genomes are around the same size and are similar overall (Globally aligned).
# Introducing integer match, mismatch, gap penalties for each iteration
# If their sizes majorly differ, massive similarity score hit
# Backtracking starts at bottom right of matrix always and matrix initialization get increasingly negative
# nw = LinearNeedlemanWunschAligner("ABxxxCDE", "ABCDE", 5, -2, -3)
# nw.execute()

# iii. Simple Fixed Smith-Waterman (Local Sequence Aligner)
# Variant of Needleman-Wunsch, finding local alignment between 2 genomic sequences.
# With SW, genmoes can be size independent as subsequence can start anywhere from the reference.
# Same match, mismatch, gap penalty scoring as NW
# Backtracking starts at whichever element gives the highest similarity score (Multiple locations possible)
# Fixed - Gap scoring is the same constant value throughout
# sw = LinearSmithWatermanAligner("AATCG", "AACG", 2, -2, -4)
# sw.execute()