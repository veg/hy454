
import multiprocessing as mp

from copy import deepcopy
from math import ceil, log
from re import compile as re_compile

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna
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
    remove_unknown = re_compile(r'[^a-zA-Z]')
    strip_front = re_compile(r'^[nN]+')
    strip_rear = re_compile(r'[nN]+$')

    for record in seqrecords:
        seq = str(record.seq)
        seq = remove_unknown.sub('', seq)
        seq = strip_front.sub('', seq)
        seq = strip_rear.sub('', seq)

        record.seq = Seq(seq, generic_dna)

    return


CUSTOM, FIRST, LONGEST = range(3)
def determine_refseq(seqrecords, mode):
    if mode not in xrange(3):
        raise ValueError("mode must be one of CUSTOM, FIRST, or LONGEST")

    if mode is CUSTOM:
        seq = raw_input("Input the entire reference sequence without newlines :: ")
        refseq = SeqRecord(Seq(seq, generic_dna),
                id="ref", name="reference",
                description="Custom reference sequence")
    elif mode is FIRST:
        refseq = seqrecords[0]
        seqrecords = seqrecords[1:]
    elif mode is LONGEST:
        idx, refseq = max(enumerate(seqrecords), key=lambda r: len(r.seq))
        seqrecords = [r for i, r in enumerate(seqrecords) if i != idx]

    return refseq, seqrecords


def _cdnaln_wrkr(refseq, seqs):
    worker = CodonAligner()
    refseqstr = str(refseq)
    forward, scores = worker.align(refseqstr, [str(s) for s in seqs])
    revcom, _scores = worker.align(refseqstr, [str(s.reverse_complement()) for s in seqs])
    assert(len(forward) == len(scores) == len(revcom) == len(_scores))
    return [forward[i] if scores[i] >= _scores[i] else revcom[i] for i in xrange(len(forward))]


def align_to_refseq(refseq, seqrecords):
    num_cpus = mp.cpu_count()
    pool = mp.Pool(num_cpus)

    seqs_per_proc = int(ceil(float(len(seqrecords)) / num_cpus))

    numseqs = len(seqrecords)
    refseqstr = str(refseq.seq)

    results = [None] * num_cpus

    for i in xrange(num_cpus):
        l = i * seqs_per_proc
        u = min(numseqs, l + seqs_per_proc)
        seqs = [s.seq for s in seqrecords[l:u]]
        results = pool.apply_async(_cdnaln_wrkr, (refseq, seqs))

    # deepcopy the seqrecords so that we can change their sequences later
    alignrecords = deepcopy(seqrecords)

    for i in xrange(num_cpus):
        l = i * seqs_per_proc
        u = min(numseqs, l + seqs_per_proc)
        seqstrs = results[i].get(0xFF)
        for j, k in enumerate(xrange(l, u)):
            alignrecords[k].seq = Seq(seqstrs[j], generic_dna)

    return MultipleSeqAlignment(alignrecords)
