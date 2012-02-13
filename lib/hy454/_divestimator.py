
import json

from collections import namedtuple
from math import ceil, log
from os.path import abspath, exists, join, split
from sys import stderr
from textwrap import dedent

from Bio.Seq import translate
from Bio.Alphabet import generic_nucleotide

from hypy import HyphyInterface, HyphyMap


__all__ = ['DiversityEstimator', 'Thresholds']


Thresholds = namedtuple('Thresholds', ['alnscore', 'minprop'])


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
        # assume MSA has been trimmed to A:B
        # strip all sequences that don't have minimum scores
        msa = [r for r, s in zip(msa, scores) if r > thresholds.alnscore]

        # bin unique variants
        amino_bins = {}
        codon_bins = {}
        for r in msa:
            codons = str(r.seq)
            aminos = translate(codons)
            # create key if non-existant
            if aminos not in amino_bins:
                amino_bins[aminos] = 0
            if codons not in codon_bins:
                codon_bins[codons] = 0
            # increment key
            amino_bins[aminos] += 1
            codon_bins[codons] += 1

        # only keep those that have the minimum copy count
        mincopy = thresholds.minprop * len(msa)
        retained = [r for r, c in amino_bins if c > mincopy]

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

        _divest_npopts = {};
        _divest_npopts[ 0 ] = "%(input_file)s";
        _divest_npopts[ 1 ] = "HKY85";
        _divest_npopts[ 2 ] = "Global";
        _divest_npopts[ 3 ] = "%(tree_file)s";
        ExecuteAFile( _divest_stdlib + "AnalyzeNucProtData.bf", _divest_npopts );
        ''') % {
            'input_file': 'blahblah',
            'tree_file': 'blahblah',
            'bootraw_file': 'blahblah',
            'bootstrap_ps_file': 'blahblah'
        })

        # self.map()

        pass
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
