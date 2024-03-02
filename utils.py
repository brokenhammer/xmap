# Useful utilities when reading files
from enum import Enum
import numpy as np

try:
    import re
    from numpy import size
except:
    print("Regular expression module 're' not available")
    raise

def file_tokens(fp):
    """ A generator to split a file into tokens
    """
    toklist = []
    while True:
        line = fp.readline()
        if not line: break
        toklist = line.split()
        for tok in toklist:
            yield tok

def file_numbers(fp):
    """Generator to get numbers from a text file"""
    toklist = []
    while True:
        line = fp.readline()
        if not line: break
        # Match numbers in the line using regular expression
        pattern = r'[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?'
        toklist = re.findall(pattern, line)
        for tok in toklist:
            yield tok

def writef(fp, val):
    """ Write values to a file in Fortran style
    
    """
    if size(val) == 1:
        fp.write("%16.9E" % val)
        return
    
    line = 0
    for v in val:
        fp.write("%16.9E" % v)
        line += 1
        if line == 5:
            fp.write("\n");
    if line != 5:
        fp.write("\n");

def np_to_list(data):
    new_data = {}
    for key,v in data.items():
        if isinstance(v, np.ndarray):
            new_data[key] = v.tolist()

        else:
            new_data[key] = v

    return new_data

class FigType(Enum):
    save = 'save'
    show = 'show'
    none = 'none'
    web = 'web'

    def __str__(self):
        return self.value