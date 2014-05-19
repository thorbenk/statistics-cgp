from matplotlib import pyplot as plot
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D

import glob
import csv
import numpy
import os

#linestyles = ['_', '-', '--', ':']
linestyles = ['-']
markers = []
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ':
            markers.append(m)
    except TypeError:
        pass
colors = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.5019607843137255, 0.5019607843137255, 0.0), (0.7529411764705882, 0.7529411764705882, 0.7529411764705882), (1.0, 0.4117647058823529, 0.7058823529411765), (0.4, 0.803921568627451, 0.6666666666666666), (0.6470588235294118, 0.16470588235294117, 0.16470588235294117), (0.0, 0.0, 0.5019607843137255), (1.0, 0.6470588235294118, 0.0), (0.6784313725490196, 1.0, 0.1843137254901961), (0.5019607843137255, 0.0, 0.5019607843137255), (0.9411764705882353, 0.9019607843137255, 0.5490196078431373)]

def mkPlot(column, ylabel, outfile, yclip=None, only=None):
    plot.clf()
    plot.figure()
    
    print "make plot", outfile
    fnames = sorted(glob.glob("stats/small_supervoxels/stat_*.txt"))
    N = 0
    for fname in fnames:
        if "NO_COMP" in fname or "ZLIB_NONE" in fname:
            continue
        
        if only is not None and only not in fname:
            continue
        
        f = open(fname, 'r')
        lines = f.readlines()
        
        r = csv.DictReader(lines, delimiter=' ')
        x = []
        y = []
        for row in r:
            x.append( float(row["sizeBytesUncompressed"]) )
            y.append( float(row[column]) )
    
        y = numpy.asarray(y);
        x = numpy.asarray(x);
        assert len(x) == len(y)
        assert x.ndim == y.ndim == 1
        
        argSort = numpy.argsort(x)
        assert argSort.ndim == 1
        assert len(argSort) == len(x)
        
        y = y[x.argsort()]
        x.sort()
        
        X = []
        Y = []
        YERR = []
    
        last = x[0]
        lastIndex = 0;
        for i, e in enumerate(x):
            if last != e:
                X.append(last)
                vals = y[lastIndex:i]
                Y.append(vals.mean())
                YERR.append(vals.std())
                lastIndex = i 
                last = e
                
        if yclip:
            plot.ylim(*yclip)
        
        X = numpy.asarray(X)
        Y = numpy.asarray(Y)
        YERR = numpy.asarray(YERR)
    
        plot.gca().set_yscale("log", nonposy='clip')
        #plot.gca().set_xscale("log", nonposx='clip')
        
        ls = numpy.random.choice(linestyles)
        c  = colors[ N % len(colors) ]
        marker = numpy.random.choice(markers)
        
        #plot.errorbar(X/1024**2.0, Y,  YERR, label=fname[len("stat_"):-4], linestyle=ls, color=c, marker=marker, markersize=3, markeredgecolor=c, markerfacecolor=c)
       
        plot.plot(X/1024**2.0, Y,  label=os.path.basename(fname)[len("stat_"):-4], linestyle=ls, color=c, marker=marker, markersize=3, markeredgecolor=c, markerfacecolor=c)
        
        plot.xlabel("uncompressed (MB)")
        plot.ylabel(ylabel)
        
        N += 1

    fontP = FontProperties()
    fontP.set_size('xx-small')
    plot.legend(prop=fontP, loc='best',borderaxespad=0.,ncol=5,mode='expand')
    plot.savefig(outfile) 

mkPlot("compessionRatio", "compression ratio", "compression_ratio_%s.pdf" % "BEST", yclip=(0,0.1), only="BEST")
mkPlot("msPerMB_compress", "ms/MB compression", "msPerMB_compress_%s.pdf" % "BEST", only="BEST")
mkPlot("msPerMB_uncompress", "ms/MB uncompression", "msPerMB_uncompress_%s.pdf" % "BEST", only="BEST")

mkPlot("compessionRatio", "compression ratio", "compression_ratio_%s.pdf" % "FAST", yclip=(0,0.4), only="FAST")
mkPlot("msPerMB_compress", "ms/MB compression", "msPerMB_compress_%s.pdf" % "FAST", only="FAST")
mkPlot("msPerMB_uncompress", "ms/MB uncompression", "msPerMB_uncompress_%s.pdf" % "FAST", only="FAST")
