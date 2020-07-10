import time
from random import randint


def compute(pFunction, pParameters, pTries):
    exception_string = ''
    for i in range(pTries):
        try:
            pFunction(pParameters)
            return
        except Exception as exp:
            time.sleep(1 + randint(0, 5))
            exception_string = str(exp)
    raise Exception('Tries {}, but still got an exception: {}'.format(pTries, exception_string))
    return
