from simpleSmithWaterman import SmithWatermanAligner

sw = SmithWatermanAligner("GAATTCAGT", "GGATCGA", 5, -2, -3)
sw.execute()
