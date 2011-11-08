
from copy import deepcopy
from math import ceil, log
from multiprocessing import cpu_count, current_process
from operator import itemgetter
from re import compile as re_compile, I as re_I
from sys import exc_info, exit as sys_exit
from types import ListType

from fakemp import create_pool

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

    if mode is CUSTOM:
        seq = raw_input("Input the entire reference sequence without newlines :: ")
        refseq = SeqRecord(Seq(seq, generic_nucleotide),
                id="ref", name="reference",
                description="Custom reference sequence")
    elif mode is FIRST:
        refseq = seqrecords[0]
        seqrecords = seqrecords[1:]
    elif mode is LONGEST:
        idx, refseq = max(enumerate(seqrecords), key=lambda r: len(r.seq))
        seqrecords = [r for i, r in enumerate(seqrecords) if i != idx]

    return refseq, seqrecords


def _cdnaln_wrkr(refseq, seqs, checkframes=xrange(3), quiet=True):
    try:
        worker = CodonAligner()
        refseqstr = str(refseq)
        bestseqs = [None] * len(seqs)
        bestscores = [0.] * len(seqs)
        for frame in checkframes:
            # zipped is a list of tuples of tuples :: [((f, fs), (r, rs)), ...]
            # worker.align return a tuple of lists, which are zipped together to get tuples of (seq, score)
            # and these are zipped to their revcom+score so that we can easily take a max later 
            zipped = zip(
                zip(worker.align(refseqstr, [str(s[frame:]) for s in seqs], quiet)),
                zip(worker.align(refseqstr, [str(s[frame:].reverse_complement()) for s in seqs], quiet))
            )
            for i in xrange(len(seqs)):
                seq, score = max(zipped[i], key=itemgetter(1))
                if score > bestscore[i]:
                    bestseqs[i] = seq
                    bestscore[i] = score
        return bestseqs
    except KeyboardInterrupt:
        return KeyboardInterrupt
    except:
        return exc_info()[1]


def align_to_refseq(refseq, seqrecords, checkframes=xrange(3)):
    num_cpus = cpu_count()

    seqs_per_proc = int(ceil(float(len(seqrecords)) / num_cpus))

    numseqs = len(seqrecords)

    try:
        results = [None] * num_cpus
        do_parts = xrange(num_cpus)
        attempts = 3
        for _ in xrange(attempts):
            pool = create_pool(refseq.seq)

            for i in do_parts:
                l = i * seqs_per_proc
                u = min(numseqs, l + seqs_per_proc)
                seqs = [s.seq for s in seqrecords[l:u]]
                results[i] = pool.apply_async(_cdnaln_wrkr, (refseq.seq, seqs, checkframes))

            pool.close()
            pool.join()

            for i in do_parts:
                results[i] = results[i].get(0xFFFF)

            if any([isinstance(r, KeyboardInterrupt) for r in results]):
                raise KeyboardInterrupt
            elif all([isinstance(r, ListType) for r in results]):
                break
            else:
                do_parts = [i for i, r in enumerate(results) if not isinstance(r, ListType)]

        excs = [e for e in results if isinstance(e, Exception)]
        if len(excs):
            raise excs[0]

        if not all([isinstance(r, ListType) for r in results]):
            raise RuntimeError("Random and unknown weirdness happened while trying to farm out work to child processes")

    except KeyboardInterrupt, e:
        if pool is not None:
            pool.terminate()
            pool.join()
        if current_process().daemon:
            return e
        else:
            print 'caught ^C (keyboard interrupt), exiting ...'
            sys_exit(-1)

    # deepcopy the seqrecords so that we can change their sequences later
    alignrecords = deepcopy(seqrecords)

    for i in xrange(num_cpus):
        l = i * seqs_per_proc
        u = min(numseqs, l + seqs_per_proc)
        for j, k in enumerate(xrange(l, u)):
            alignrecords[k].seq = Seq(results[i][j], generic_nucleotide)

    return MultipleSeqAlignment(alignrecords)
