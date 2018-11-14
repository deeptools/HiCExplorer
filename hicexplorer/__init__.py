import logging
logging.basicConfig(level=logging.DEBUG)
import warnings 
import sys

if not sys.warnoptions:
    warnings.simplefilter("ignore")