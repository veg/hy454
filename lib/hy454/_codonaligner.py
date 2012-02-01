
import json

from math import ceil, log
from os.path import abspath, exists, join, split
from sys import stderr

from Bio.Alphabet import generic_nucleotide

from hypy import HyphyMap


__all__ = ['CodonAligner']


class CodonAligner(HyphyMap):

    def __init__(self, batchfile=None, retvar=None):
        if batchfile is None:
            batchfile = join(
                    split(abspath(__file__))[0],
                    'hyphy', 'codonaligner.bf'
            )

        if retvar is None:
            retvar = "_cdnaln_outstr"

        if not exists(batchfile):
            raise ValueError("Invalid batchfile `%s', it doesn't exist!" % batchfile)
        super(CodonAligner, self).__init__(batchfile, retvar)

    def __call__(self, refseq, seqs, revcomp=False, quiet=True):
        return CodonAligner.align(self, refseq, seqs, revcomp, quiet)

    def align(self, refseq, seqs, revcomp=False, quiet=True):
        # if we have no sequences, abort early to prevent later errors
        if not len(seqs):
            return [], []

        # uppercase the refseq to deal with bugs in HyPhy's codon aligner
        refseq = refseq.upper()

        # pad the reference to the nearest codon,
        # otherwise the hyphy codon alignment algo barfs 
        if len(refseq) > 3:
            pad = 3 - (len(refseq) % 3)
            pad = 0 if pad == 3 else pad
            refseq += '-' * pad
            scoremod = float(len(refseq)) / (len(refseq) - pad)
        else:
            pad = 0
            scoremod = 0.

        numseqs = len(seqs)
        # if the # nodes exceeds the number of seqs, we just need numseqs jobs
        numnodes = min(numseqs, self.nodes)
        seqs_per_node = max(1, numseqs // numnodes)
        remainder = numseqs % numnodes

        arg1 = 'Yes' if revcomp else 'No'

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
            argslist.append( [arg1, refseq, len(node_seqs)] + node_seqs )

        retstrs = self.map(argslist, quiet)

        seqscores = []
        for retstr in retstrs:
            seqscores.extend(json.loads(retstr))

        newseqstrs, scores = zip(*seqscores)

        return list(newseqstrs), list(scores)
