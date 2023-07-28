import numpy as np
import sys

ics_dynspecs = [np.load(fname) for fname in sys.argv[2:]]
sum_ics = sum(ics_dynspecs)
np.save(sys.argv[1], sum_ics)
