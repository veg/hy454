#
# idepi :: (IDentify EPItope) python libraries containing some useful machine
# learning interfaces for regression and discrete analysis (including
# cross-validation, grid-search, and maximum-relevance/mRMR feature selection)
# and utilities to help identify neutralizing antibody epitopes via machine
# learning.
#
# Copyright (C) 2011 N Lance Hepler <nlhepler@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide, generic_protein


__all__ = ['OrfList']


class OrfList(object):

    def __init__(self, seq):
        # accumulate all the potential start codon indices
        seq = seq.upper()
        start_idxs = []
        for codon in ('ATG', 'AUG'):
            prefix = 0
            while True:
                idx = seq.find(codon)
                if idx >= 0:
                    start_idxs.append(idx+prefix)
                    idx += 3
                    seq = seq[idx:]
                    prefix += idx
                else:
                    break

        start_idxs = sorted(start_idxs)

        orfs = []

        # find only those that terminate.. and truncate to its stop codon
        dnaseq = Seq(seq, generic_nucleotide)

        for idx in start_idxs:
            orfseq = dnaseq[idx]
            ter_idx = orfseq.translate().find('*')
            # if we don't find a termination codon, go ahead and add the full thing to the list
            if ter_idx < 0:
                pass
            # cut out small peptides (those <10 amino acids long)
            elif ter_idx < 10:
                continue
            else:
                orfseq = orfseq[:(3 * ter_idx)]
            orfs.append(orfseq)

        self.__orfs = sorted(orfs, key=lambda x: len(x), reverse=True)

    def __getitem__(self, key):
        return self.__orfs[key]

    def __len__(self):
        return len(self.__orfs)

    def __contains__(self, key):
        return True if key >= 0 and key < len(self.__orfs) else False
