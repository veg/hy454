
from __future__ import division, print_function

import json

from os.path import abspath, exists, join, split

from BioExt import (BLOSUM62, DNAExpIdScoreMatrix,
    DNAScoreMatrix, ProteinScoreMatrix)

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

        self.__codon = codon

        super(Aligner, self).__init__(batchfile, retvar)

    @property
    def codon(self):
        return self.__codon

    def __call__(self, refseq, seqs, score_matrix=None, revcomp=False, expected_identity=0., keep_insertions=False, quiet=True):
        return Aligner.align(self, refseq, seqs, score_matrix, revcomp, expected_identity, keep_insertions, quiet)

    def align(self, refseq, seqs, score_matrix=None, revcomp=False, expected_identity=0., keep_insertions=False, quiet=True):
        # if we have no sequences, abort early to prevent later errors
        if not len(seqs):
            return [], []

        if score_matrix is None:
            if self.codon:
                score_matrix = BLOSUM62.load()
            else:
                score_matrix = DNAExpIdScoreMatrix(
                    0.8 if expected_identity == 0. else expected_identity,
                    { 'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25 }
                )

        if self.codon and not isinstance(score_matrix, ProteinScoreMatrix):
            raise ValueError('score_matrix incompatible with codon alignment')
        elif not self.codon and not isinstance(score_matrix, DNAScoreMatrix):
            raise ValueError('score_matrix incompatible with dna alignment')

        smdef = { ('_cdnaln_letters', '_cdnaln_scorematrix'): score_matrix }

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

        retstrs = self.map(argslist, globalvars=smdef, quiet=quiet)

        seqscores = []
        for retstr in retstrs:
            seqscores.extend(json.loads(retstr))

        newrefstrs, newseqstrs, scores, overlaps, identities = zip(*seqscores)

        return list(newrefstrs), list(newseqstrs), list(scores), list(overlaps), list(identities)
