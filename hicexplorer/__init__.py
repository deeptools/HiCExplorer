import logging
# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('cooler').setLevel(logging.ERROR)
logging.getLogger('hicmatrix').setLevel(logging.DEBUG)
logging.getLogger('numexpr').setLevel(logging.ERROR)


import warnings
import sys

if not sys.warnoptions:
    warnings.simplefilter("ignore")

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
