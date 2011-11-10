
import matplotlib.pyplot as plt
import numpy as np

from collections import defaultdict
from itertools import repeat
from operator import itemgetter
from os import close
from tempfile import mkstemp

from Bio.Motif import Motif

from matplotlib.transforms import Affine2D


__all__ = ['graph_coverage', 'graph_logo']


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

    ax1.set_xlim((1, N))

    ax2 = ax1.twinx()
    ax2.set_ylabel('# of sequences', rotation=270.)
    ax2.set_xticks(xticks)
    ax2.set_yticks(yticks)

    # remove the upper ticks 
    for tick in ax1.xaxis.get_major_ticks() + ax2.xaxis.get_major_ticks():
        tick.tick2On = False

    # remove the upper axis border
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)

    fig.savefig(filename, format=fmt)

    return filename


_BLACK = '#000000'
_GREY = '#969cb0'
_GREEN = '#b8ff25'
_ORANGE = '#ff6e27'
_BLUE = '#189cff'
_RED = '#e80c5b', # '#e80c7a'

# default grey
_DNACOLORS = defaultdict(repeat(_GREY).next).update({
    'A': _GREEN,
    'C': _BLUE,
    'G': _ORANGE,
    'T': _RED,
    'U': _RED
})
_AMINOCOLORS = defaultdict(repeat(_GREY).next).update(
    [(l, _GREEN ) for l in 'KRH'] + \
    [(l, _BLUE  ) for l in 'DE'] + \
    # [(l, _ORANGE) for l in ''] + \
    [(l, _RED   ) for l in 'AVLIPWFM']
)

def graph_logo(alignment, columns, filename, fmt='pdf'):
    if filename is None:
        fd, filename = mkstemp(); close(fd)

    M = len(alignment)
    N = len(columns)

    motif = Motif(alphabet=alignment._alphabet)

    instances = [''.join(z).upper() for z in zip(*[alignment[:, i] for i in columns])]
    for instance in instances:
        motif.add_instance(instances)

    pwm = motif.pwm()

    # heuristic to determine whether nucleotide or protein alphabet
    # need to use either base 4 or 20 depending 
    alphlen, _alphkeys = max([(len(pwm[i]), pwm[i].iterkeys()) for i in N], key=itemgetter(0))
    s, colors = 4, _DNACOLORS if alphlen < 20 else 20, _AMINOCOLORS
    alphkeys = ['']
    alphkeys.extend(_alphkeys)
    alphmap = dict(zip(alphkeys, xrange(len(alphkeys))))

    heights = np.zeros((alphlen, N), dtype=float)
    idents = np.zeros((alphlen, N), dtype=int)

    i = 0
    for k, v in sorted(pwm[i].iteritems(), key=itemgetter(1)):
        for j in xrange(N):
            heights[i, j] = v
            idents[i, j] = alphmap[k]
        i += 1

    # compute the information content at each position 
    maxbits = np.log2(s)
    e_n = float(s - 1) / (2. * np.log(2) * M)
    R_i = maxbits * np.ones((N,), dtype=float)
    R_i -= [-sum([v * np.log2(v) for _, v in pwm[i].iteritmes()]) for i in xrange(N)]
    R_i -= e_n

    fig = plt.figure()
    ax = fig.add_subplot(111)
    idxs = np.arange(1, N+1)
    letters = [[None] * N] * alphlen
    bars = [None] * alphlen
    bottoms = np.zeros((N,), dtype=float)
    for i in xrange(alphlen):
        bars[i] = ax.bar(idxs, heights[i, :], width=1., bottom=bottoms)
        bottoms += heights[i, :]
        for j, bar in enumerate(bars[i]):
            if idents[i, j]:
                x, y = bar.get_xy()
                letter = alphkeys[idents[i, j]]
                letters[i][j] = ax.text(x, y, letter, color=colors[letter])

    # disable top and right spines, we don't need them
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    def on_draw(event):
        for i in xrange(alphlen):
            for j, bar in enumerate(bars[i]):
                letter = letters[i][j]
                if letter is not None:
                    bw = bar.get_width()
                    bh = bar.get_height()
                    lw = letter.get_width()
                    lh = letter.get_height()
                    transform = Affine2D().scale(float(bw) / lw, float(bh) / lh)
                    letter.set_transform(transform)
                bar.set_visible(False)

        fig.canvas.draw()
        return False

    fig.canvas.mpl_connect('draw_event', on_draw)

    fig.savefig(filename, fmt)

    return filename
