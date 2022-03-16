import logging
logging.basicConfig(level=logging.DEBUG)
# logging.basicConfig(level=logging.INFO)
logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('cooler').setLevel(logging.ERROR)
logging.getLogger('hicmatrix').setLevel(logging.ERROR)
# logging.getLogger('hicmatrix').setLevel(logging.DEBUG)
logging.getLogger('hicmatrix').setLevel(logging.ERROR)

logging.getLogger('numba').setLevel(logging.ERROR)

import os
os.environ["OMP_NUM_THREADS"] = "10" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "10" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "10" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "10" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "10" # export NUMEXPR_NUM_THREADS=6
os.environ["NUMEXPR_MAX_THREADS"] = "10" # export NUMEXPR_NUM_THREADS=6



import warnings
import sys

if not sys.warnoptions:
    warnings.simplefilter("ignore")

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

import multiprocessing as mp
mp.set_start_method('fork')
