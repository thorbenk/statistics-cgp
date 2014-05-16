from matplotlib import pyplot as plot
from matplotlib.font_manager import FontProperties

import glob
import csv
import numpy

fnames = sorted(glob.glob("stat_*.txt"))
for fname in fnames:
    if "NO_COMP" in fname or "ZLIB_NONE" in fname:
        continue
    print "*** file =", fname,"***"
    
    f = open(fname, 'r')
    lines = f.readlines()
    
    r = csv.DictReader(lines, delimiter=' ')
    x = []
    compressionRatio = []
    for row in r:
        print row
        x.append( float(row["sizeBytesUncompressed"]) )
        compressionRatio.append( float(row["compessionRatio"]) )
   
    compressionRatio = numpy.asarray(compressionRatio);
    x = numpy.asarray(x);
    assert len(x) == len(compressionRatio)
    assert x.ndim == compressionRatio.ndim == 1
    
    argSort = numpy.argsort(x)
    assert argSort.ndim == 1
    assert len(argSort) == len(x)
    
    print type(compressionRatio), compressionRatio.shape, compressionRatio.dtype
    print type(x),  x.shape, x.dtype
    
    #from IPython import embed; embed()
    
    y = compressionRatio[x.argsort()]
    x.sort()
    
    plot.plot(x/1024**2.0,compressionRatio.copy(), "-", label=fname[len("stat_"):-4])
    plot.xlabel("uncompressed (MB)")
    plot.ylabel("compression ratio")

fontP = FontProperties()
fontP.set_size('small')
plot.legend(prop=fontP)

plot.show()



