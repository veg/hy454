
from math import ceil, log
from os.path import abspath, exists, join, split
from sys import stderr

from Bio.Alphabet import generic_nucleotide

from hypy import HyphyInterface


__all__ = ['CodonAligner']


class CodonAligner(HyphyInterface):

    def __init__(self, batchfile=None):
        if batchfile is None:
            batchfile = join(
                    split(abspath(__file__))[0],
                    'hyphy', 'codonaligner.bf'
            )
        if not exists(batchfile):
            raise ValueError("Invalid batchfile `%s', it doesn't exist!" % batchfile)
        # use only 1 cpu
        super(CodonAligner, self).__init__(batchfile, 1)

    def align(self, refseq, seqs, quiet=True, return_by_position=False):
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

        # compute next power of two size string for the output
        outlen = (len(refseq) + 1) * len(seqs)
        outlen = ceil(log(outlen, 2))
        outlen = int(pow(2, outlen))

        self.queuestralloc('_cdnaln_outstr', outlen)
        self.queuevar('_cdnaln_refseq', refseq)
        self.queuevar('_cdnaln_seqs', seqs)
        if return_by_position:
        	self.queuevar('_cdnaln_returnByPosition', 1.0)
       	

        self.runqueue()

        if not quiet:
            if self.stdout != '':
                print(self.stdout, file=stderr)
            if self.warnings != '':
                print(self.warnings, file=stderr)

        if self.stderr != '':
            raise RuntimeError(self.stderr)

        newseqstrs = self.getvar('seqs', HyphyInterface.STRING).split(',')
        
        scores = self.getvar('scores', HyphyInterface.MATRIX)
        if pad:
            newseqstrs = [s[:-pad] for s in newseqstrs]
        if scoremod:
            scores = [s * scoremod for s in scores]
        return newseqstrs, scores
