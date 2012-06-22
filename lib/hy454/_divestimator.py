
from __future__ import division, print_function

import json

from collections import Counter, namedtuple
from math import ceil, log
from os import close, remove
from os.path import abspath, exists, join, split
from sys import stderr
from tempfile import mkstemp
from textwrap import dedent

from Bio.Seq import translate
from Bio.Alphabet import generic_nucleotide

from hppy import HyphyInterface, HyphyMap


__all__ = ['DiversityEstimator', 'Thresholds']


class Thresholds(object):

    def __init__(self, alnscore, mincopy=None, minprop=None):
        if mincopy == None and minprop == None:
            raise ValueError('one of either mincopy or minprop must be provided')
        self.alnscore = alnscore
        self.mincopy = 0 if mincopy is None else mincopy
        self.minprop = 0 if minprop is None else minprop


class DiversityEstimator(HyphyMap):

    def __init__(self, batchfile=None, retvar=None):
        if batchfile is None:
            batchfile = join(
                    split(abspath(__file__))[0],
                    'hyphy', 'divestimator.bf'
            )

        if retvar is None:
            retvar = "_divest_outstr"

        if not exists(batchfile):
            raise ValueError("Invalid batchfile `%s', it doesn't exist!" % batchfile)
        super(DiversityEstimator, self).__init__(batchfile, retvar)

    def __call__(self, refseq, msa, scores, thresholds, quiet=True):
        return DiversityEstimator.estimate(self, refseq, msa, scores, quiet)

    def estimate(self, refseq, msa, scores, thresholds, quiet=True):
        if not isinstance(thresholds, Thresholds):
            raise TypeError("unsupported operand type(s) for estimate: '%s'" % str(type(thresholds)))

        # assume MSA has been trimmed to A:B
        # bin unique variants
        # strip all sequences that don't have minimum scores
        nuc_bins = Counter(str(r.seq) for r, s in zip(msa, scores) if s > thresholds.alnscore)
        amino_bins = Counter()
        for cdn, c in nuc_bins.items():
            amino_bins[translate(cdn)] += c

        numvariants = sum(nuc_bins.values())

        # only keep those that have the minimum copy count
        mincopy = max(thresholds.mincopy, thresholds.minprop * numvariants)
        retained = [r for r, c in amino_bins if c > mincopy]

        inputfile, treefile, bootrawfile, bspsfile = [None] * 4
        try:
            fd, inputfile = mkstemp(); close(fd)
            fd, treefile = mkstemp(); close(fd)
            fd, bootrawfile = mkstemp(); close(fd)
            fd, bspsfile = mkstemp(); close(fd)

            with open(inputfile, 'w') as fh:
                print('\n'.join('>%d\n%s' % (i, seq) for i, seq in enumerate(retained)), file=fh)

            iface = HyphyInterface()
            iface.runqueue(execstr=dedent('''\
            _divest_stdlib = HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR;
    
            _divest_njopts = {};
            _divest_njopts[ 0 ] = "Distance formulae";
            _divest_njopts[ 1 ] = "Nucleotide/Protein";
            _divest_njopts[ 2 ] = "%(input_file)s";
            _divest_njotps[ 3 ] = "Force Zero";
            _divest_njopts[ 4 ] = "TN93";
            _divest_njopts[ 5 ] = "y";
            _divest_njopts[ 6 ] = "%(tree_file)s";
            ExecuteAFile( _divest_stdlib + "NeighborJoining.bf", _divest_njopts );
    
            _divest_bsopts = {};
            _divest_bsopts[ 0 ] = "100";
            _divest_bsopts[ 1 ] = "%(bootraw_file)s";
            _divest_bsopts[ 2 ] = "Proportions";
            _divest_bsopts[ 3 ] = "70";
            _divest_bsopts[ 4 ] = "%(bootstrap_ps_file)s";
            ExecuteAFile( _divest_stdlib + "post_npbs.bf", _divest_bsopts );
    
            _divest_anpdopts = {};
            _divest_anpdopts[ 0 ] = "%(input_file)s";
            _divest_anpdopts[ 1 ] = "HKY85";
            _divest_anpdopts[ 2 ] = "Global";
            _divest_anpdopts[ 3 ] = "%(tree_file)s";
            ExecuteAFile( _divest_stdlib + "AnalyzeNucProtData.bf", _divest_anpdopts );
            ''') % {
                'input_file': inputfile,
                'tree_file': treefile,
                'bootraw_file': bootrawfile,
                'bootstrap_ps_file': bspsfile
            })

            # self.map()
        finally:
            if exists(inputfile):
                remove(inputfile)
            if exists(treefile):
                remove(treefile)
            if exists(bootrawfile):
                remove(bootrawfile)
            if exists(bspsfile):
                remove(bspsfile)

        # define window A, B
        # align sequences to ref, truncating to [A, B)
        # create bins for each obs in window by aa, and codons
        # the total # of unique variants == len(binner)
        # filter variants by count or % threshold (* len(sequences))
        # do NJ with "Distance formulae", "Nucleotide/protein", inFile, "Force Zero", "TN93", "y", inFile + ".tree"
        #   @ HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "NeighborJoining.bf"
        # do bootstrap with "100", inFile + ".bootraw", "Proportions", "70", inFile + ".njboot.ps"
        #   @ HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "post_npbs.bf"
        # do NucProt analysis with inFile, "HKY85", "Global", inFile + ".tree"
        #   @ HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "AnalyzeNucProtData.bf" 
        # for each replicate (Map part):
        #   do NJ like above with 2: inFile + ".sim" and 6: inFile + ".sim.tree", and "nj" in EAB
        #   do NucProt analysis like above with 0: inFile + ".sim" and 3: inFile + ".sim.tree", and "sim" in EAB
        #   do simD[replicate] = maxPathInATree( "sim.givenTree" ); 
