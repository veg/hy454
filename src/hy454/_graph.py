
import matplotlib.pyplot as plt
import numpy as np

from collections import defaultdict
from itertools import repeat
from operator import itemgetter
from os import close
from tempfile import mkstemp

from Bio.Alphabet import Gapped, HasStopCodon, _verify_alphabet
from Bio.Alphabet.IUPAC import \
        ambiguous_dna, ambiguous_rna, extended_protein, \
        unambiguous_dna, unambiguous_rna
from Bio.Motif import Motif
from Bio.Seq import Seq

from matplotlib.transforms import Affine2D


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

    ax1.fill_between(np.arange(1, N+1), height, edgecolor=_GREY, facecolor=_GREY, alpha=0.5)
    # ax1.bar(np.arange(1, N+1), height, width=1., color='black', edgecolor='black')
    ax1.set_xlabel('Reference sequence position')
    ax1.set_ylabel('Coverage')
    # we don't need to set the xticks here because we do it for ax2 
    # ax1.set_xticks(xticks)
    ax1.set_yticks(np.arange(0.2, 1.1, 0.2))

    ax1.set_xlim((0.5, N+1.5))

    ax2 = ax1.twinx()
    ax2.set_ylabel('# of sequences', rotation=270.)
    ax2.set_xticks(xticks+0.5)
    ax2.set_xticklabels([str(int(t)) for t in xticks])
    ax2.set_yticks(yticks)

    # remove the upper ticks 
    for tick in ax1.xaxis.get_major_ticks() + ax2.xaxis.get_major_ticks():
        tick.tick2On = False

    # remove the upper axis border
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)

    fig.savefig(filename, format=fmt)

    return filename


_DNA_ALPHABET = Gapped(ambiguous_dna)
_RNA_ALPHABET = Gapped(ambiguous_rna)
_AMINO_ALPHABET = HasStopCodon(Gapped(extended_protein))

def _fix_ambigs(pwm, alphabet):
    # remove ambigs by distributing their probability uniformly 
    if alphabet == _DNA_ALPHABET or alphabet == _RNA_ALPHABET:
        T = 'T' if alphabet == _DNA_ALPHABET else 'U'
        mapper = {
            'M': ('A', 'C'),
            'R': ('A', 'G'),
            'W': ('A',  T ),
            'S': ('C', 'G'),
            'Y': ('C',  T ),
            'K': ('G',  T ),
            'V': ('A', 'C', 'G'),
            'H': ('A', 'C',  T ),
            'D': ('A', 'G',  T ),
            'B': ('C', 'G',  T ),
            'N': ('A', 'C', 'G', T )
        }
        for i in xrange(len(pwm)):
            for k, unambig in mapper.iteritems():
                if k in pwm[i]:
                    # redistribute the probability uniformly
                    C = pwm[i][k] / len(unambig)
                    for l in unambig:
                        pwm[i][l] += C
                    # remove the key
                    del pwm[i][k]
    # remove gaps
    for i in xrange(len(pwm)):
        if _GAP in pwm[i]:
            # uniform redistribution to everybody else, eg p(letter | observed)
            C = pwm[i][_GAP] / (len(pwm[i])-1)
            for l in pwm[i].iterkeys():
                pwm[i][l] += C
            # remove the gap
            del pwm[i][_GAP]
    return pwm


_BLACK = '#000000'
_GREY = '#969cb0'
_GREEN = '#b8ff25'
_ORANGE = '#ff6e27'
_BLUE = '#189cff'
_RED = '#e80c5b' # '#e80c7a'

# default grey
_DNA_COLORS = defaultdict(repeat(_GREY).next, {
    'A': _GREEN,
    'C': _BLUE,
    'G': _ORANGE,
    'T': _RED,
    'U': _RED
})
_AMINO_COLORS = defaultdict(repeat(_GREY).next,
    [(l, _GREEN ) for l in 'KRH'] + \
    [(l, _BLUE  ) for l in 'DE'] + \
    # [(l, _ORANGE) for l in ''] + \
    [(l, _RED   ) for l in 'AVLIPWFM']
)

_LOGO_UPDATED = False

def graph_logo(alignment, columns, filename, fmt='pdf'):
    global _LOGO_UPDATED

    if filename is None:
        fd, filename = mkstemp(); close(fd)

    M = len(alignment)
    N = len(columns)

    alph = None
    for _alph in (_DNA_ALPHABET, _RNA_ALPHABET, _AMINO_ALPHABET):
        for r in alignment:
            r.seq.alphabet = _alph
        if all([_verify_alphabet(r.seq) for r in alignment]):
            alph = _alph
            break
    if alph is None:
        raise RuntimeError("sequences with indeterminable alphabet provided")

    motif = Motif(alphabet=alph)

    instances = [''.join(z).upper() for z in zip(*[alignment[:, i] for i in columns])]
    for instance in instances:
        motif.add_instance(Seq(instance, alph))

    pwm = _fix_ambigs(motif.pwm(), alph)

    # heuristic to determine whether nucleotide or protein alphabet
    # need to use either base 4 or 20 depending 
    alphlen, _alphkeys = max([(len(pwm[i]), pwm[i].iterkeys()) for i in xrange(N)], key=itemgetter(0))
    s, colors = (4, _DNA_COLORS) if alphlen < 20 else (20, _AMINO_COLORS)
    alphkeys = ['']
    alphkeys.extend(_alphkeys)
    alphmap = dict(zip(alphkeys, xrange(len(alphkeys))))

    # compute the information content at each position 
    maxbits = np.log2(s)
    e_n = float(s - 1) / (2. * np.log(2) * M)
    R = maxbits * np.ones((N,), dtype=float)
    R -= [-sum([v * np.log2(v) for _, v in pwm[i].iteritems()]) for i in xrange(N)]
    R -= e_n

    heights = np.zeros((alphlen, N), dtype=float)
    idents = np.zeros((alphlen, N), dtype=int)

    for j in xrange(N):
        i = 0
        for k, v in sorted(pwm[i].iteritems(), key=itemgetter(1)):
            heights[i, j] = R[j] * v
            idents[i, j] = alphmap[k]
            i += 1

    _LOGO_UPDATED = False
    fig = plt.figure()
    ax = fig.add_subplot(111)
    idxs = np.arange(1, N+1)
    barletters = [[None] * N] * alphlen
    bottoms = np.zeros((N,), dtype=float)
    for i in xrange(alphlen):
        bars = ax.bar(idxs, heights[i, :], width=1., bottom=bottoms)
        bottoms += heights[i, :]
        for j, bar in enumerate(bars):
            if idents[i, j]:
                x, y = bar.get_xy()
                l = alphkeys[idents[i, j]]
                barletters[i][j] = bar, ax.text(x, y, l, color=colors[l])

    ax.set_ylim((0, maxbits))

    # remove the top and right ticks 
    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.tick2On = False

    # disable top and right spines, we don't need them
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks(np.arange(1, N+1, dtype=int))

    def on_draw(event):
        global _LOGO_UPDATED
        if not _LOGO_UPDATED:
            for i in xrange(alphlen):
                for j in xrange(N):
                    bar, letter = barletters[i][j]
                    x, y, bw, bh = bar.get_window_extent().bounds
                    _, _, lw, lh = letter.get_window_extent().bounds
                    print x, y, bw, bh, lw, lh
                    tr = Affine2D().scale(100, 100)
                    tr = tr.translate(x, y)
                    letter.set_transform(tr)
                    bar.set_visible(False)
            _LOGO_UPDATED = True
            fig.canvas.draw()
        return False

    fig.canvas.mpl_connect('draw_event', on_draw)

    fig.show()

    fig.savefig(filename, format=fmt)

    return filename
