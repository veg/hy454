
from __future__ import division, print_function

import matplotlib as mpl
mpl.use('pdf')

import matplotlib.pyplot as plt
import numpy as np

from collections import defaultdict, Counter
from itertools import repeat
from operator import itemgetter
from os import close
from os.path import dirname, join
from re import compile as re_compile
from tempfile import mkstemp

from Bio.Alphabet import Gapped, HasStopCodon, _verify_alphabet
from Bio.Alphabet.IUPAC import (ambiguous_dna, ambiguous_rna,
        extended_protein, unambiguous_dna, unambiguous_rna)
from Bio.Motif import Motif
from Bio.Seq import Seq

from matplotlib.font_manager import (FontProperties,
        createFontList, findSystemFonts, fontManager)
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from matplotlib.transforms import Affine2D

from BioExt import _GAP, _STOP

from ._basefont import Basefont


__all__ = [
    'COVERAGE',
    'MAJORITY',
    'graph_coverage_majority',
    'graph_logo'
]


_HY454_FONT_PATHS = [join(dirname(__file__), 'data', 'fonts', 'ttf')]


# update the fontManager to handle the Roboto font installed with hy454
# fontManager.ttffiles.extend(findSystemFonts(_HY454_FONT_PATHS))
# fontManager.ttflist = createFontList(fontManager.ttffiles)
_ROBOTO_REGULAR = FontProperties(fname=join(_HY454_FONT_PATHS[0], 'Roboto-Regular.ttf'))


def _adjust_spines_outward(ax, spines, points):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', points))
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')


def _max_nonzero_min(values, default=0):
    vmax = max(values)
    vmin = vmax
    for v in values:
        if v > 0 and v < vmin:
            vmin = v
    vmin = default if vmin == vmax else vmin
    return vmax, vmin


COVERAGE, MAJORITY = 1, 2

def graph_coverage_majority(
    alignment,
    mode,
    filename=None,
    dpi=None, figsize=None, format='pdf', transparent=True
):
    if not mode:
        mode = COVERAGE | MAJORITY

    if mode not in (COVERAGE, MAJORITY, COVERAGE | MAJORITY):
        raise ValueError('mode requires COVERAGE, MAJORITY, or COVERAGE | MAJORITY')

    if filename is None:
        fd, filename = mkstemp(); close(fd)

    if figsize is None:
        # actually results in a graph that is almost 4" high
        # and 6.4" wide -- almost the golden ratio
        figsize = (6, 6)

    n0 = 0
    M = len(alignment)
    N = alignment.get_alignment_length()

    # if we don't match the whole reference,
    # cut off the head, tail of gaps
    leading = re_compile(r'^-*')
    trailing = re_compile(r'-*$')
    l = -N
    n = -N
    for r in alignment:
        seq = str(r.seq)
        # *ing.search should always match
        v1 = len(leading.search(seq).group(0))
        v2 = len(trailing.search(seq).group(0))
        if v1 < abs(l):
            l = v1
        if v2 < abs(n):
            n = v2
    # if we're greater than 0, remove leading, trailing gaps
    if l > 0:
        n0 += l
    if n > 0:
        N -= n

    frac = 1. / M
    heights = np.zeros((N - n0,), dtype=float)

    xs = np.arange(n0 + 1, N + 1)

    # some heuristic to prevent too many ticks,
    # is complicated by the start and end dynamically business
    # 0.477122 jumps an order 10 at 30% of the next order 10
    specialK = 0.477122
    xdiv = 10 ** int(np.log10(N) - specialK) or 1
    xsep = xdiv * int((N - n0) / 5 / xdiv + 1)
    xstart = (n0 // xsep) * xsep
    xend = N + xsep
    xticks = np.arange(xstart, xend, xsep)
    if len(xticks) > 1:
        fix = False
        delta0 = xticks[1] - n0
        deltaN = N - xticks[-2]
        min_delta = xsep / 4
        if delta0 != 0 and delta0 < min_delta:
            xstart += xsep
            fix = True
        if deltaN != 0 and deltaN < min_delta:
            xend -= xsep
            fix = True
        if fix:
            xticks = np.arange(xstart, xend, xsep)
    xticks[0] = n0 + 1
    xticks[-1] = N

    ydiv = 10 ** int(np.log10(M) - specialK) or 1
    ysep = ydiv * int(M / 5 / ydiv + 1)
    yticks = np.arange(ysep, M+ysep, ysep)
    if len(yticks) > 1:
        if M - yticks[-2] < ysep / 4:
            yticks = np.arange(ysep, M, ysep)
    yticks[-1] = M

    # mpl.rcParams['font.family'] = 'sans-serif'
    # mpl.rcParams['font.sans-serif'] = 'Roboto'
    # mpl.rcParams['font.weight'] = '400'

    fig = plt.figure(figsize=figsize, dpi=dpi)

    # golden rectangle!
    rect = 0.2, 0.2, 1, 0.618
    ax1 = fig.add_axes(rect)

    # move the axis spines off the data
    # _adjust_spines_outward(ax1, ('bottom', 'left', 'right'), 18)

    majorities = np.zeros((N - n0,), dtype=float)
    if mode & MAJORITY:
        for i, col in enumerate(range(n0, N)):
            # count the number of occurrences of each unique character
            # then grab the maximum count observed and multiply it by
            # frac to get the majorities
            counts = Counter(alignment[:, col].upper())
            if _GAP in counts:
                del counts[_GAP]
            # grab the count of the most common variant
            m = 0
            if len(counts):
                m = counts.most_common(1)[0][1]
            # compute the normalized majority in the majority-only case,
            # taking care to avoid division-by-zero errors
            if mode == MAJORITY:
                s = sum(counts.values())
                s = 1 if s == 0 else s
                majorities[i] = m / s
            else:
                majorities[i] = m * frac
        ax1.plot(
            xs, majorities,
            color=_LRED, linewidth=1., zorder=-1, alpha=0.8 if mode & COVERAGE else 1.0
        )

    if mode & COVERAGE:
        for i, col in enumerate(range(n0, N)):
            heights[i] = sum(frac for p in alignment[:, col] if p != _GAP)

    # labels
    ax1.set_xlabel('Reference sequence position', fontproperties=_ROBOTO_REGULAR)

    extra_artists = []

    if mode == (COVERAGE | MAJORITY):
        ax1.plot(
            xs, heights,
            color=_LBLUE, linewidth=1., zorder=-2
        )
        # create a proxy artist for legend, PolyCollections don't work (heights)
        p1 = Line2D([0, 1], [0, 1], color=_LBLUE, linewidth=1.)
        # create a proxy artist for legend, [Lines2D] don't work (majorities)
        p2 = Line2D([0, 1], [0, 1], color=_LRED, linewidth=1.)
        leg = ax1.legend(
            [p1, p2], ['Coverage', 'Majority'],
            bbox_to_anchor=(0.5, -0.15), loc=9, ncol=2,
            prop=_ROBOTO_REGULAR, borderpad=0.
        )
        leg.legendPatch.set_alpha(0.)
        extra_artists.append(leg)
        # ax1.set_ylabel('Coverage - majority', fontproperties=_ROBOTO_REGULAR)
    elif mode == COVERAGE:
        ax1.fill_between(
            xs, heights, majorities,
            edgecolor=_LBLUE, facecolor=_LBLUE, linewidth=1., zorder=-1
        )
        ax1.set_ylabel('Coverage', fontproperties=_ROBOTO_REGULAR)
    else:
        ax1.set_ylabel('Proportion', fontproperties=_ROBOTO_REGULAR)

    def format_percent(x, pos=None):
        return '%1.0f%%' % (100 * x)

    ax1.yaxis.set_major_formatter(FuncFormatter(format_percent))
    ax1.set_yticks(np.arange(0.2, 1.1, 0.2))
    ax1.set_xlim((n0 + 1, N))
    ax1.set_ylim((0, 1))

    major_ticks = ax1.xaxis.get_major_ticks(len(xticks))

    # disable the first and last tick on the x-axis,
    # they're redundant and ugly (esp if minH > 0)
    major_ticks[0].tick1On = False
    major_ticks[-1].tick1On = False

    # if we're only doing coverage, include the number
    # of sequences, but with majority in the mix
    # it doesn't make sense
    if mode == MAJORITY:
        ax1.set_xticks(xticks)
        # if we're not showing the # of sequences on the left,
        # remove the ticks and spine
        major_ticks += ax1.yaxis.get_major_ticks()
        ax1.spines['right'].set_visible(False)
        # alter axes to show max, min value
        maxM, minM = _max_nonzero_min(majorities)
        ax1.spines['left'].set_bounds(minM, maxM)
        # set font properties
        ticklabels = (
            ax1.xaxis.get_ticklabels() +
            ax1.yaxis.get_ticklabels()
        )
    else:
        ax2 = ax1.twinx()
        # move the axis spines off the data
        # _adjust_spines_outward(ax2, ('bottom', 'left', 'right'), 18)
        ax2.set_ylabel('No. of sequences', rotation=270., fontproperties=_ROBOTO_REGULAR)
        ax2.set_xticks(xticks)
        ax2.set_yticks(yticks)
        # set transparent here otherwise ax2 doesn't exist
        if transparent:
            ax2.patch.set_alpha(0.)
        # get the major ticks so we can disable them later
        major_ticks += ax2.xaxis.get_major_ticks()
        # disable the top spines, like we do later
        ax2.spines['top'].set_visible(False)
        # use the axes spines to show the maximum value
        maxH, minH = _max_nonzero_min(heights)
        ax1.spines['left'].set_bounds(minH, maxH)
        ax1.spines['right'].set_bounds(minH, maxH)
        # set font properties
        ticklabels = (
            ax1.xaxis.get_ticklabels() +
            ax1.yaxis.get_ticklabels() +
            ax2.xaxis.get_ticklabels() +
            ax2.yaxis.get_ticklabels()
        )

    for label in ticklabels:
        label.set_fontproperties(_ROBOTO_REGULAR)

    if transparent:
        fig.patch.set_alpha(0.)
        ax1.patch.set_alpha(0.)

    # remove the upper ticks
    for tick in major_ticks:
        tick.tick2On = False

    # remove the upper axis border
    ax1.spines['top'].set_visible(False)

    fig.savefig(
        filename,
        format=format,
        transparent=transparent,
        bbox_extra_artists=extra_artists,
        bbox_inches='tight',
        pad_inches=0.25
    )

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


def graph_logo(
    alignment,
    columns,
    filename=None,
    dpi=None, edgecolor='k', figsize=None, format='pdf', labels=None, linewidth=0., transparent=True
):
    if filename is None:
        fd, filename = mkstemp(); close(fd)

    if figsize is None:
        figsize = (3, 3)

    if labels is None:
        labels = ['%d' % (idx + 1) for idx in columns]

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

    # set laplace = True to include the backgrounds
    pwm = _fix_ambigs(motif.pwm(laplace=False), alph)

    # heuristic to determine whether nucleotide or protein alphabet
    # need to use either base 4 or 20 depending
    alphlen, _alphkeys = max([(len(pwm[i]), pwm[i].keys()) for i in range(N)], key=itemgetter(0))
    s, colors = (4, _DNA_COLORS) if alphlen < 20 else (20, _AMINO_COLORS)
    alphkeys = ['']
    alphkeys.extend(_alphkeys)
    alphmap = dict(zip(alphkeys, range(len(alphkeys))))

    # compute the information content at each position
    maxbits = np.log2(s)
    e_n = (s - 1) / (2. * np.log(2) * M)
    R = maxbits * np.ones((N,), dtype=float)
    R -= [-sum([v * np.log2(v) for _, v in pwm[i].items() if v > 0.]) for i in range(N)]
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

    # mpl.rcParams['font.family'] = 'sans-serif'
    # mpl.rcParams['font.sans-serif'] = 'Roboto'
    # mpl.rcParams['font.weight'] = '400'

    fig = plt.figure(figsize=figsize, dpi=dpi)

    # make each column a vertical golden rect
    rect = 0.2, 0.2, 0.382 * N, 0.618
    ax = fig.add_axes(rect)

    ax.set_ylabel('bits', fontproperties=_ROBOTO_REGULAR)

    if figsize is None:
        fig.set_figwidth(N)

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

    # set font properties
    for label in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
        label.set_fontproperties(_ROBOTO_REGULAR)

    # disable top and right spines, we don't need them
    ax.spines['bottom'].set_visible(False)
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
                glyph = font[l]
                ax.add_patch(glyph)
                glyph.set_transform(bar.get_transform())
                bar.set_visible(False)
                glyph.set_edgecolor(edgecolor)
                glyph.set_facecolor(colors[l])
                glyph.set_linewidth(linewidth)
                glyph.set_zorder(-1)

    # set the remaining spine to show the maximum value
    ax.spines['left'].set_bounds(0., max(bottoms))

    fig.savefig(filename, format=format, transparent=transparent, bbox_inches='tight', pad_inches=0.25)

    return filename
