import time
from random import randint
import traceback


def compute(pFunction, pParameters, pTries):
    exception_string = ''
    pTries = 5
    for i in range(pTries):
        try:
            pFunction(pParameters)
            return
        except Exception:
            time.sleep(1 + randint(0, 5))
            exception_string = traceback.format_exc()
    raise Exception('Tries {}, but still got an exception: {}'.format(pTries, exception_string))
    return
