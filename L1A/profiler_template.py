# Run a profiler on select code

import cProfile
import pstats
import pdb

cProfile.run('import camal_l1a_v0',filename='profiler_results')
st = pstats.Stats('profiler_results')
st.sort_stats('calls')
st.print_stats(0.03)

pdb.set_trace()
