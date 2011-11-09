
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
    height = np.zeros((N,), dtype=float)
    sep = int(N / 5 / 25) * 25

    for i in xrange(N):
        height[i] = sum([frac for s in alignment if s.seq[i] != _GAP])

    plt.bar(np.arange(N), height, width=1., color='black', edgecolor='black')
    plt.xlabel('Reference Sequence Position')
    plt.ylabel('Coverage')
    plt.xticks(np.arange(0, N+1, sep))
    plt.yticks(np.arange(0, 1.01, 0.2))
    plt.savefig(filename, format=fmt)

    return filename
