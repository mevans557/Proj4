import numpy as np
import scipy as sp
import sys

sys.settrace

A = sp.sparse.random_array((10,10))
print(A)
