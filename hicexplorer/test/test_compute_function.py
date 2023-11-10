import time
from random import randint
import traceback

def compute(pFunction, pParameters, pTries):
    exception_string = ''
    pTries = 1
    for i in range(pTries):
        try:
            pFunction(pParameters)
            return
        except Exception as exp:
            time.sleep(1 + randint(0, 5))
            # detailed_traceback = 
            exception_string = traceback.format_exc()
    raise Exception('Tries {}, but still got an exception: {}'.format(pTries, exception_string))
    return
