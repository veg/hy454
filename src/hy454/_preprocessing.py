
from copy import deepcopy
from math import ceil, log
from multiprocessing import cpu_count, current_process
from operator import itemgetter
from re import compile as re_compile, I as re_I
from sys import exc_info, exit as sys_exit, float_info
from types import ListType

from fakemp import farmout, farmworker

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_nucleotide
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from _codonaligner import CodonAligner


__all__ = [
    'preprocess_seqrecords',
    'CUSTOM', 'FIRST', 'LONGEST',
    'determine_refseq',
    'align_to_refseq'
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
    if mode not in xrange(3):
        raise ValueError("mode must be one of CUSTOM, FIRST, or LONGEST")

    if mode == CUSTOM:
        seq = raw_input("Input the entire reference sequence without newlines :: ")
        refseq = SeqRecord(Seq(seq, generic_nucleotide),
                id="ref", name="reference",
                description="Custom reference sequence")
    elif mode == FIRST:
        refseq = seqrecords.pop(0)
    elif mode == LONGEST:
        idx, refseq = max(enumerate(seqrecords), key=lambda r: len(r[1].seq))
        seqrecords.pop(idx)

    return refseq, seqrecords


def _codonaligner(refseq, seqs, quiet=True):
    if not len(seqs):
        return []
    worker = CodonAligner()
    refseqstr = str(refseq)
    # zipped is a list of tuples of tuples :: [((f, fs), (r, rs)), ...]
    # worker.align return a tuple of lists, which are zipped together to get tuples of (seq, score)
    # and these are zipped to their revcom+score so that we can easily take a max later 
    zipped = zip(
        zip(*worker.align(refseqstr, [str(s) for s in seqs], quiet)),
        zip(*worker.align(refseqstr, [str(s.reverse_complement()) for s in seqs], quiet))
    )
    # itemgetter(1) is the score, [0] is the seq
    return [max(z, key=itemgetter(1))[0] for z in zipped]


def _worker(refseq, seqs, quiet=True):
    return farmworker(_codonaligner, refseq, seqs, quiet)


def align_to_refseq(refseq, seqrecords):
    num_cpus = cpu_count()

    seqs_per_proc = int(ceil(float(len(seqrecords)) / num_cpus))

    numseqs = len(seqrecords)

    results = farmout(
        num_cpus,
        lambda i: (refseq.seq, [s.seq for s in seqrecords[(i*seqs_per_proc):min(numseqs, (i+1)*seqs_per_proc)]]),
        _worker,
        lambda r: isinstance(r, ListType),
        attempts=3
    )

    # deepcopy the seqrecords so that we can change their sequences later
    alignrecords = deepcopy(seqrecords)

    for i in xrange(num_cpus):
        l = i * seqs_per_proc
        u = min(numseqs, l + seqs_per_proc)
        for j, k in enumerate(xrange(l, u)):
            alignrecords[k].seq = Seq(results[i][j], generic_nucleotide)

    return MultipleSeqAlignment(alignrecords)
