
import matplotlib.pyplot as plt
import numpy as np

from os import close
from tempfile import mkstemp

from Bio.Motif import Motif


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

    height = [sum([frac for p in alignment[:, i] if p != _GAP]) for i in xrange(N)]

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

    # remove the upper ticks 
    for tick in ax1.xaxis.get_major_ticks() + ax2.xaxis.get_major_ticks():
        tick.tick2On = False

    # remove the upper axis border
    ax1.spines['top'].cla() # set_linewidth(0)
    ax2.spines['top'].cla() # set_linewidth(0)

    fig.savefig(filename, format=fmt)

    return filename


def graph_logo(alignment, columns, filename, fmt='pdf'):
    if filename is None:
        fd, filename = mkstemp(); close(fd)

    M = len(alignment)
    N = len(columns)

    motif = Motif(alphabet=alignment._alphabet)

    instances = [''.join(z) for z in zip(*[alignment[:, i] for i in columns])]
    for instance in instances:
        motif.add_instance(instances)

    pwm = motif.pwm()

    # heuristic to determine whether nucleotide or protein alphabet
    # need to use either base 4 or 20 depending 
    s = 4 if len(pwm[0]) < 20 else 20

    # compute the information content at each position 
    maxbits = np.log2(s)
    e_n = float(s - 1) / (2. * np.log(2) * M)
    R_i = maxbits * np.ones((N,), dtype=float)
    R_i -= [-sum([v * np.log2(v) for _, v in pwm[i].iteritmes()]) for i in xrange(N)]
    R_i -= e_n
