

from copy import deepcopy
from json import dump as json_dump
from math import ceil, log
from multiprocessing import cpu_count, current_process
from operator import itemgetter
from re import compile as re_compile, I as re_I
from sys import exc_info, exit as sys_exit, float_info

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_nucleotide
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BioExt import enumerate_by_codon

from fakemp import farmout, farmworker

from ._codonaligner import CodonAligner
from ._graph import _GAP


__all__ = [
    'preprocess_seqrecords',
    'CUSTOM', 'FIRST', 'LONGEST',
    'determine_refseq',
    'align_to_refseq',
    'positional_write'
]


def preprocess_seqrecords(seqrecords):
    remove_unknown = re_compile(r'[^ACGTUWSMKRYBDHVN]', re_I)
    strip_front = re_compile(r'^[N]+', re_I)
    strip_rear = re_compile(r'[N]+$', re_I)

    for record in seqrecords:
        seq = str(record.seq)
        seq = remove_unknown.sub('', seq)
        seq = strip_front.sub('', seq)
        seq = strip_rear.sub('', seq)

        record.seq = Seq(seq, generic_nucleotide)

    return


CUSTOM, FIRST, LONGEST = range(3)
def determine_refseq(seqrecords, mode):
    if mode not in range(3):
        raise ValueError("mode must be one of CUSTOM, FIRST, or LONGEST")

    if mode == CUSTOM:
        seq = input("Input the entire reference sequence without newlines :: ")
        refseq = SeqRecord(Seq(seq, generic_nucleotide),
                id="ref", name="reference",
                description="Custom reference sequence")
    elif mode == FIRST:
        refseq = seqrecords.pop(0)
    elif mode == LONGEST:
        idx, refseq = max(enumerate(seqrecords), key=lambda r: len(r[1].seq))
        seqrecords.pop(idx)

    return refseq, seqrecords


def align_to_refseq(refseq, seqrecords, revcomp=True, quiet=False):

    aligned, scores = CodonAligner()(str(refseq.seq), [str(s.seq) for s in seqrecords], revcomp, quiet)

    # deepcopy the seqrecords so that we can change their sequences later
    alignrecords = deepcopy(seqrecords)

    for i, aln in enumerate(aligned):
        alignrecords[i].seq = Seq(aln, generic_nucleotide)

    return MultipleSeqAlignment(alignrecords)


def positional_write(msa, fh):
    datastruct = {}
    for seq in msa:
        seqdata = []
        for pos, cdn in enumerate_by_codon(seq):
            if cdn != (_GAP * 3):
                seqdata.append((pos, cdn))
        if len(seqdata):
            datastruct[seq.id] = seqdata

    json_dump(datastruct, fh)
    fh.write('\n')
