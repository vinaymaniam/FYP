import cProfile
import re
import pstats

from runUFRESH import runUFRESH

import cProfile, pstats, io
pr = cProfile.Profile()
pr.enable()

runUFRESH()

pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats(10)
print(s.getvalue())
