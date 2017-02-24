import numpy as np
import sys 
import ast 
import pyximport; pyximport.install()
import duplicateRemover 
import ast
import warnings


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--sum', dest='METHOD', action='store_const',const="sum", default="max",  help='use sum of distances between two sides of the read (default: use max distance)')
parser.add_argument("--mindist", type=int, default = 3,  help = 'Reads less than this distance (in bp) away are considered duplicates. Default: 3')
parser.add_argument("--sep", type=str, default=r"\v", help=r"Separator (\t, \v, etc. characters are supported, pass them in quotes) ")
parser.add_argument("--out", type=str, default="", help="File to dump duplicates, default stdout")
parser.add_argument("--dupfile", type=str, default="", help="File to dump duplicates, if desired. By default they are discarded")

parser.add_argument("--c1", type=int, default = 0,  help = 'Chrom 1 column; default 0')
parser.add_argument("--c2", type=int, default = 3,  help = 'Chrom 2 column; default 3')
parser.add_argument("--p1", type=int, default = 1,  help = 'Position 1 column; default 1')
parser.add_argument("--p2", type=int, default = 4,  help = 'Position 2 column; default 4')
parser.add_argument("--s1", type=int, default = 2,  help = 'Strand 1 column; default 2')
parser.add_argument("--s2", type=int, default = 5,  help = 'Strand 2 column; default 5')

args = parser.parse_args()

sep = ast.literal_eval('"""' + args.sep + '"""')

METHOD = args.METHOD
MINDIST = args.mindist

c1ind = args.c1
c2ind = args.c2
p1ind = args.p1
p2ind = args.p2
s1ind = args.s1 
s2ind = args.s2

outfile = args.out
if outfile:
    outstream = open(outfile,'w')
else:
    outstream = sys.stdout

otherfile = args.dupfile
if otherfile:
    dupfile = open(otherfile,'w')

maxind = max(c1ind, c2ind, p1ind, p2ind, s1ind, s2ind)

# you don't need to load more than 10k lines at a time b/c you get out of the CPU cache, so this parameter is not adjustable
maxLen = 10000 

myclass = duplicateRemover.processiveDuplicateRemoval(METHOD, MINDIST, returnData=False)


def fetchadd(key, mydict):
    key = key.strip()
    if key not in mydict:
        mydict[key] = len(mydict)
    return mydict[key]

def ar(mylist, val):
    dtype = {8:np.int8, 16:np.int16, 32:np.int32}[val]
    return np.array(mylist, dtype=dtype)
    

c1 = []; c2 = []; p1 = []; p2 = []; s1 = []; s2 = []
lines = []
chromDict = {}
strandDict = {}


while True: 
    line = sys.stdin.readline()    
    stripline = line.strip()
    if stripline.startswith("#"):
        outstream.write(line)
        continue
    
    if line:
        if not stripline: 
            warnings.warn("Empty line detected not at the end of the file")
            continue            
        lines.append(line)
        words = line.split(sep)
        if len(words) <= maxind:
            raise ValueError("Error parsing line {0}: expected {1} words, got {2}".format(line, maxind, len(words)))
            
        c1.append(fetchadd(words[c1ind], chromDict))
        c2.append(fetchadd(words[c2ind], chromDict))
        p1.append(int(words[p1ind]))
        p2.append(int(words[p2ind]))
        s1.append(fetchadd(words[s1ind], strandDict))
        s2.append(fetchadd(words[s2ind], strandDict))
    if (not line) or (len(c1) == maxLen):
        res = myclass.addChunk(ar(c1,8), ar(c2,8), ar(p1,32), ar(p2, 32), ar(s1, 8), ar(s2,8))
        if not line:
            res = np.concatenate([res, myclass.finish()])
        for newline, remove in zip(lines[:len(res)], res):
            if not remove:
                outstream.write(newline)  
            else:
                if otherfile:
                    dupfile.write(newline)
                
        c1 = []; c2 = []; p1 = []; p2 = []; s1 = []; s2 = []
        lines = lines[len(res):]
        if not line:
            if(len(lines) != 0):                
                raise ValueError("{0} lines left in the buffer, should be none; something went terribly wrong".format(len(lines)))
            break 
