
import matplotlib.pyplot as plt
import numpy as np

from os import close
from tempfile import mkstemp


__all__ = ['graph_coverage']


_GAP = '-'

def graph_coverage(refseq, alignment, filename=None, fmt='pdf'):
    if filename is None:
        fd, filename = mkstemp(); close(fd)

    N = len(refseq)
    frac = 1. / len(alignment)
    height = np.zeros((N,), dtype=int)

    for i in xrange(N):
        height[i] = round(sum([frac for s in alignment if s[i] != _GAP]))

    plt.bar(np.arange(N), height, width=1., color='black', edgecolor='black')
    plt.xlabel('Reference Sequence Position')
    plt.ylabel('Coverage (%)')
    plt.xticks(np.arange(0, N+1, 20))
    plt.yticks(np.arange(0, 101, 20))
    plt.savefig(filename, format=fmt)

    return filename
