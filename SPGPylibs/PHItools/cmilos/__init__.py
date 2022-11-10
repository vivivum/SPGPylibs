try:
    from .pymilos import pmilos
except:
    print("unable to import pmilos in __init__.py in .PHItools.cmilos (this is o.k.)")

# This is important for pdoc to ignore other files, like setup.py
__all__ = ['pmilos']
