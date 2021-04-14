import functools
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
