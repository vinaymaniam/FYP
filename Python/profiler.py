import cProfile
import re
import pstats

from full_training import runFullTraining

import cProfile, pstats, io

def profile():
    pr = cProfile.Profile()
    pr.enable()

    runFullTraining(0,1,0)
    pr.disable()
    s = io.StringIO()
    sortby = 'tottime'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(10)
    print(s.getvalue())
    return
profile()
