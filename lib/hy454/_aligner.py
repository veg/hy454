
from __future__ import division, print_function

import json

from math import ceil, log
from os.path import abspath, exists, join, split
from sys import stderr

from Bio.Alphabet import generic_nucleotide

from hppy import HyphyMap


__all__ = ['Aligner']


class Aligner(HyphyMap):

    def __init__(self, codon=True, batchfile=None, retvar=None):
        if batchfile is None:
            batchfile = join(
                    split(abspath(__file__))[0],
                    'hyphy', 'codonaligner.bf' if codon else 'aligner.bf'
            )

        if retvar is None:
            retvar = '_cdnaln_outstr' if codon else '_aln_outstr'

        if not exists(batchfile):
            raise ValueError("Invalid batchfile `%s', it doesn't exist!" % batchfile)
        super(Aligner, self).__init__(batchfile, retvar)

    def __call__(self, refseq, seqs, revcomp=False, expected_identity=0., keep_insertions=False, quiet=True):
        return Aligner.align(self, refseq, seqs, revcomp, expected_identity, keep_insertions, quiet)

    def align(self, refseq, seqs, revcomp=False, expected_identity=0., keep_insertions=False, quiet=True):
        # if we have no sequences, abort early to prevent later errors
        if not len(seqs):
            return [], []

        # uppercase the refseq to deal with bugs in HyPhy's aligner
        refseq = refseq.upper()

        numseqs = len(seqs)
        # if the # nodes exceeds the number of seqs, we just need numseqs jobs
        numnodes = min(numseqs, self.nodes)
        seqs_per_node = max(1, numseqs // numnodes)
        remainder = numseqs % numnodes

        arg1 = 'Yes' if revcomp else 'No'
        arg2 = 'Yes' if keep_insertions else 'No'

        argslist = []
        lwr, upr = 0, 0

        for i in range(numnodes):
            # since our traversal is stateful, keep these cursors
            # around. During the first remainder iterations,
            # add an extra seq to the list of seqs, afterwards
            # proceed as normal
            lwr = upr
            if i < remainder:
                upr = min(numseqs, lwr + seqs_per_node + 1)
            else:
                upr = min(numseqs, lwr + seqs_per_node)
            node_seqs = [s.upper() for s in seqs[lwr:upr]]
            argslist.append( [arg1, arg2, refseq, expected_identity, len(node_seqs)] + node_seqs )

        retstrs = self.map(argslist, quiet)

        seqscores = []
        for retstr in retstrs:
            seqscores.extend(json.loads(retstr))

        newseqstrs, scores, overlaps, identities = zip(*seqscores)

        return list(newseqstrs), list(scores), list(overlaps), list(identities)
