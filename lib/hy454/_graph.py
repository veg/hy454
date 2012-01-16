
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from collections import defaultdict
from itertools import repeat
from operator import itemgetter
from os import close
from os.path import dirname, join
from tempfile import mkstemp

from Bio.Alphabet import Gapped, HasStopCodon, _verify_alphabet
from Bio.Alphabet.IUPAC import (ambiguous_dna, ambiguous_rna,
        extended_protein, unambiguous_dna, unambiguous_rna)
from Bio.Motif import Motif
from Bio.Seq import Seq

from matplotlib.font_manager import (createFontList,
        findSystemFonts, fontManager)
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from matplotlib.transforms import Affine2D

from ._basefont import Basefont


__all__ = ['graph_coverage', 'graph_logo']


_GAP = '-'
_STOP = '*'


_HY454_FONT_PATHS = [join(dirname(__file__), 'data', 'fonts', 'ttf')]


# update the fontManager to handle the Roboto font installed with hy454
fontManager.ttffiles.extend(findSystemFonts(_HY454_FONT_PATHS))
fontManager.ttflist = createFontList(fontManager.ttffiles)


def graph_coverage(alignment, filename=None, format='pdf', transparent=True):
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

    height = [sum([frac for p in alignment[:, i] if p != _GAP]) for i in range(N)]

    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'Roboto'

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.fill_between(np.arange(1, N+1), height, edgecolor=_LBLUE, facecolor=_LBLUE, linewidth=0., zorder=-1)
    ax1.set_xlabel('Reference sequence position')
    ax1.set_ylabel('Coverage')
    # we don't need to set the xticks here because we do it for ax2 
    # ax1.set_xticks(xticks)
    ax1.set_yticks(np.arange(0.2, 1.1, 0.2))

    ax1.set_xlim((1., N))

    ax2 = ax1.twinx()
    ax2.set_ylabel('# of sequences', rotation=270.)
    ax2.set_xticks(xticks)
    ax2.set_xticklabels([str(int(t)) for t in xticks])
    ax2.set_yticks(yticks)

    if transparent:
        fig.patch.set_alpha(0.)
        ax1.patch.set_alpha(0.)
        ax2.patch.set_alpha(0.)

    # remove the upper ticks 
    for tick in ax1.xaxis.get_major_ticks() + ax2.xaxis.get_major_ticks():
        tick.tick2On = False

    # remove the upper axis border
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)

    fig.savefig(filename, format=format, transparent=transparent)

    return filename


_DNA_ALPHABET = Gapped(ambiguous_dna)
_RNA_ALPHABET = Gapped(ambiguous_rna)
_AMINO_ALPHABET = HasStopCodon(Gapped(extended_protein, gap_char=_GAP), stop_symbol=_STOP)

def _fix_ambigs(pwm, alphabet):
    mapper = {}
    killchars = _GAP
    # killchars ambigs by distributing their probability uniformly 
    if alphabet == _DNA_ALPHABET or alphabet == _RNA_ALPHABET:
        T = 'T' if alphabet == _DNA_ALPHABET else 'U'
        mapper.update({
            'M': 'AC',
            'R': 'AG',
            'W': 'A' + T,
            'S': 'CG',
            'Y': 'C' + T,
            'K': 'G' + T,
            'V': 'ACG',
            'H': 'AC' + T,
            'D': 'AG' + T,
            'B': 'CG' + T,
            # 'N': 'ACG' + T
        })
        killchars += 'N'
    elif alphabet == _AMINO_ALPHABET:
        mapper.update({
            'B': 'DN',
            'J': 'IL',
            'Z': 'EQ',
            # 'X': 'ACDEFGHIKLMNPQRSTVWY',
        })
        killchars += _STOP + 'X' + 'OU'
    for i in range(len(pwm)):
        for k, unambig in mapper.items():
            if k in pwm[i]:
                # redistribute the probability uniformly
                C = pwm[i][k] / len(unambig)
                for l in unambig:
                    pwm[i][l] += C
                # killchars the key
                del pwm[i][k]
    # killchars gaps
    for i in range(len(pwm)):
        for char in killchars:
            if char in pwm[i]:
                # uniform redistribution to everybody else, eg p(letter | observed)
                C = pwm[i][char] / (len(pwm[i])-1)
                for l in pwm[i].keys():
                    pwm[i][l] += C
                # remove the gap
                del pwm[i][char]
    return pwm


_BLACK = '#000000'
_GREY = '#969cb0'
_GREEN = '#b8ff25'
_ORANGE = '#ff6e27'
_BLUE = '#189cff'
_RED = '#e80c5b' # '#e80c7a'

# android color swatch
_LGREY = '#F2F2F2'
_LBLUE = '#33B5E5'
_LPURPLE = '#AA66CC'
_LGREEN = '#99CC00'
_LORANGE = '#FFBB33'
_LRED = '#FF4444'
_DGREY = '#DDDDDD'
_DBLUE = '#0099CC'
_DPURPLE = '#9933CC'
_DGREEN = '#669900'
_DORANGE = '#FF8800'
_DRED = '#CC0000'

# default grey
_DNA_COLORS = defaultdict(repeat(_GREY).__next__, {
    'A': _LGREEN,
    'C': _LBLUE,
    'G': _LORANGE,
    'T': _LRED,
    'U': _LRED
})
_AMINO_COLORS = defaultdict(repeat(_GREY).__next__,
    [(l, _LGREEN ) for l in 'KRH'] + \
    [(l, _LBLUE  ) for l in 'DE'] + \
    # [(l, _ORANGE) for l in ''] + \
    [(l, _LRED   ) for l in 'AVLIPWFM']
)


def graph_logo(alignment, columns, filename, dpi=None, edgecolor='k', figsize=None, format='pdf', labels=None, linewidth=0., transparent=True):
    if filename is None:
        fd, filename = mkstemp(); close(fd)

    if labels is None:
        labels = ['%d' % idx for idx in columns]

    M = len(alignment)
    N = len(columns)

    alph = None
    for _alph in (_DNA_ALPHABET, _RNA_ALPHABET, _AMINO_ALPHABET):
        for r in alignment:
            r.seq.alphabet = _alph
        if all([_verify_alphabet(r.seq.upper()) for r in alignment]):
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
    alphlen, _alphkeys = max([(len(pwm[i]), pwm[i].keys()) for i in range(N)], key=itemgetter(0))
    s, colors = (4, _DNA_COLORS) if alphlen < 20 else (20, _AMINO_COLORS)
    alphkeys = ['']
    alphkeys.extend(_alphkeys)
    alphmap = dict(zip(alphkeys, range(len(alphkeys))))

    # compute the information content at each position 
    maxbits = np.log2(s)
    e_n = float(s - 1) / (2. * np.log(2) * M)
    R = maxbits * np.ones((N,), dtype=float)
    R -= [-sum([v * np.log2(v) for _, v in pwm[i].items()]) for i in range(N)]
    R -= e_n

    heights = np.zeros((alphlen, N), dtype=float)
    identities = np.zeros((alphlen, N), dtype=int)

    for j in range(N):
        i = 0
        for k, v in sorted(pwm[j].items(), key=itemgetter(1)):
            heights[i, j] = R[j] * v
            identities[i, j] = alphmap[k]
            i += 1

    font = Basefont(join(_HY454_FONT_PATHS[0], 'Roboto-Black.ttf'))

    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'Roboto'

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    if transparent:
        fig.patch.set_alpha(0.)
        ax.patch.set_alpha(0.)

    # remove the top and right ticks
    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.tick2On = False

    # remove the bottom ticks
    for tick in ax.xaxis.get_major_ticks():
        tick.tick1On = False

    # rotate the x-axis labels by 45 degrees to enhance packing
    for label in ax.xaxis.get_ticklabels():
        label.set_rotation(45)

    # disable top and right spines, we don't need them
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    def format_xlabel(x, pos=None):
        idx = np.clip(int(x)-1, 0, N-1)
        return labels[idx]

    ax.xaxis.set_major_formatter(FuncFormatter(format_xlabel))
    # avoid too much precision
    ax.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))

    # set the ticks
    ax.set_yticks(np.append(np.arange(0, maxbits, 0.5, dtype=float), maxbits))
    ax.set_xticks(np.arange(1, N+1, dtype=float) + 0.5)

    # set the axes limits here AFTER the ticks, otherwise borkage
    ax.set_xlim((1, N+1))
    ax.set_ylim((0, maxbits))

    idxs = np.arange(1, N+1)
    bottoms = np.zeros((N,), dtype=float)
    for i in range(alphlen):
        bars = ax.bar(idxs, heights[i, :], width=1., bottom=bottoms)
        bottoms += heights[i, :]
        for j, bar in enumerate(bars):
            if identities[i, j]:
                l = alphkeys[identities[i, j]]
                glyph = font.char_patch(l)
                ax.add_patch(glyph)
                glyph.set_transform(bar.get_transform())
                bar.set_visible(False)
                glyph.set_edgecolor(edgecolor)
                glyph.set_facecolor(colors[l])
                glyph.set_linewidth(linewidth)
                glyph.set_zorder(-1)

    fig.savefig(filename, format=format, transparent=transparent)

    return filename
