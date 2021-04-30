import functools

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    RESET = '\u001b[0m'

def printc(*args, color = bcolors.RESET, **kwargs):
    """My custom print() function."""
    print(u"\u001b"+f"{color}", end='\r')
    print(*args, **kwargs)
    print(u"\u001b"+f"{bcolors.RESET}", end='\r')
    return 

def countcalls(fn):
    "decorator function count function calls. Use @countcalls above def "

    @functools.wraps(fn)
    def wrapped(*args):
        wrapped.ncalls +=1
        return fn(*args)

    wrapped.ncalls = 0
    return wrapped

import time
def timeit(method):
    '''helper function to calculate executing time. Use @timeit above def'''
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            t = (te - ts) * 1000
            if t > 1e3:
                print('%r  %2.2f seconds' %
                      (method.__name__, t/1000.))
            elif t > 1e6:
                print('%r  %2.2f mim' %
                      (method.__name__, t/1000./60.))
            else:
                print('%r  %2.2f ms' %
                      (method.__name__, t))

        return result
    return timed
