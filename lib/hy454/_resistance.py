
from os.path import abspath, exists, join, split
from sys import stderr
from itertools import product

__all__ = ['DRAM','codon_mapper_max', 'codon_mapper_min', 'codon_mapper_avg']

_codonTable = {'AAA' : 'K', 'AAC' : 'N', 'AAG' : 'K', 'AAT' : 'N', 'ACA' : 'T', 'ACC' : 'T', 'ACG' : 'T', 'ACT' : 'T', 'AGA' : 'R', 'AGC' : 'S', 'AGG' : 'R', 'AGT' : 'S', 'ATA' : 'I', 'ATC' : 'I', 'ATG' : 'M', 'ATT' : 'I', 'CAA' : 'Q', 'CAC' : 'H', 'CAG' : 'Q', 'CAT' : 'H', 'CCA' : 'P', 'CCC' : 'P', 'CCG' : 'P', 'CCT' : 'P', 'CGA' : 'R', 'CGC' : 'R', 'CGG' : 'R', 'CGT' : 'R', 'CTA' : 'L', 'CTC' : 'L', 'CTG' : 'L', 'CTT' : 'L', 'GAA' : 'E', 'GAC' : 'D', 'GAG' : 'E', 'GAT' : 'D', 'GCA' : 'A', 'GCC' : 'A', 'GCG' : 'A', 'GCT' : 'A', 'GGA' : 'G', 'GGC' : 'G', 'GGG' : 'G', 'GGT' : 'G', 'GTA' : 'V', 'GTC' : 'V', 'GTG' : 'V', 'GTT' : 'V', 'TAA' : 'X', 'TAC' : 'Y', 'TAG' : 'X', 'TAT' : 'Y', 'TCA' : 'S', 'TCC' : 'S', 'TCG' : 'S', 'TCT' : 'S', 'TGA' : 'X', 'TGC' : 'C', 'TGG' : 'W', 'TGT' : 'C', 'TTA' : 'L', 'TTC' : 'F', 'TTG' : 'L', 'TTT' : 'F'}
_nucCodes   = {'A' : ['A'], 'C' : ['C'], 'G' : ['G'], 'T' : ['T'], 'M' : ['A','C'], 'R' : ['A', 'G'], 'W' : ['A', 'T'], 'S' : ['C', 'G'], 'Y' : ['C','T'], 'K' : ['G', 'T'], 'V': ['A','C','G'], 'H' : ['A','C','T'], 'D' : ['A','G','T'], 'B' : ['C', 'G', 'T']}

def codon_mapper_max(codon):
    #codon_mapped = translate_codon (codon)
    pass

def translate_codon(codon):
    if codon in _codonTable:
        return [_codonTable[codon]]
    else:
        try:
            position_translations = [_nucCodes[letter] for letter in codon]
            amino_acids = set()
            for resolution in product (position_translations[0], position_translations[1], position_translations[2]):
                amino_acids.add(_codonTable[''.join(resolution)])

            return list(amino_acids)

        except KeyError:
            return []

COMPLETE_LIST, BINARY = range(2)

def DRAM (position_list, codon_mapper=None, report_mode=None):
    if codon_mapper is None:
        codon_mapper = codon_mapper_max
    if report_mode is None:
        report_mode = COMPLETE_LIST
    pass

