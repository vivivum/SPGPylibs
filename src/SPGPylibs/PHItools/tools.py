import os, time, functools


class bcolors:
    """
    This is a simple class for colors.
    """
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    YELLOW = '\033[93m'
    WARNING = '\033[36m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    RESET = '\u001b[0m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'


def printc(*args, color=bcolors.RESET, **kwargs):
    """
    This function wraps the python ``print()`` functions adding color capabilities.

    :param args: params to pass through to ``print()`` function.
    :param color: provide the text color
    :param kwargs: **kwargs enables printc to retain all ``print()`` functionality
    :return: nothing
    """
    print(u"\u001b" + f"{color}", end='\r')
    print(*args, **kwargs)
    print(u"\u001b" + f"{bcolors.RESET}", end='\r')
    return


def countcalls(fn):
    """
    Decorator function count function calls. Use @countcalls above def.
    """

    @functools.wraps(fn)
    def wrapped(*args):
        wrapped.ncalls += 1
        return fn(*args)

    wrapped.ncalls = 0
    return wrapped


def timeit(method):
    """
    Helper function to calculate executing time. Use @timeit above def.
    """

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
                      (method.__name__, t / 1000.))
            elif t > 1e6:
                print('%r  %2.2f mim' %
                      (method.__name__, t / 1000. / 60.))
            else:
                print('%r  %2.2f ms' %
                      (method.__name__, t))

        return result

    return timed


def fix_path(path, dir='forward', verbose=False):
    """
    Corrects the path from UNIX filesystem \) -> ), etc...
    :param path: Path to fix
    :param dir: default is forward
    :param verbose: print path
    :return: corrected path
    """
    path = repr(path)
    if dir == 'backward':
        path = path.replace("\\\\", "")
        path = path.split("'")[1]
        if verbose == True:
            print('backward')
            print(path)
        return path
    elif dir == 'forward':
        path = path.replace(")", "\)")
        path = path.replace("(", "\(")
        path = path.replace(" ", "\ ")
        path = os.path.abspath(path).split("'")[1]
        if verbose == True:
            print('forward')
            print(path)
        return path


def find_div(x: int, lim: int) -> None:
    """
    Finds the divisible numbers of x
    :param x: input number
    :type x: integer
    :param lim: find divisible up to lim
    :return: None
    """

    def divisible(m, n):
        return m % n == 0

    for i in range(lim):
        if divisible(x, i + 1):
            print(i)
