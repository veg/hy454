
import matplotlib.pyplot as plt
import numpy as np

from os import close
from tempfile import mkstemp


__all__ = ['graph_coverage']


_GAP = '-'

def graph_coverage(alignment, filename=None, fmt='pdf'):
    if filename is None:
        fd, filename = mkstemp(); close(fd)

    M = len(alignment)
    N = alignment.get_alignment_length()
    frac = 1. / M
    height = np.zeros((N,), dtype=float)

    xdiv = 10 ** int(np.log10(N)-0.3) or 1
    xsep = int(float(N) / 5 / xdiv + 1) * xdiv
    xticks = np.arange(0, N+xsep, xsep)
    xticks[0] = 1
    xticks[-1] = N

    ydiv = 10 ** int(np.log10(M)-0.3) or 1
    ysep = int(float(M) / 5 / ydiv + 1) * ydiv
    yticks = np.arange(ysep, M+ysep, ysep)
    yticks[-1] = M

    for i in xrange(N):
        height[i] = sum([frac for s in alignment if s.seq[i] != _GAP])

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.bar(np.arange(1, N+1), height, width=1., color='black', edgecolor='black')
    ax1.set_xlabel('Reference sequence position')
    ax1.set_ylabel('Coverage')
    # we don't need to set the xticks here because we do it for ax2 
    # ax1.set_xticks(xticks)
    ax1.set_yticks(np.arange(0.2, 1.1, 0.2))

    ax2 = ax1.twinx()
    ax2.set_ylabel('# of sequences', rotation=270.)
    ax2.set_xticks(xticks)
    ax2.set_yticks(yticks)

    fig.savefig(filename, format=fmt)

    return filename
