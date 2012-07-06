
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BioExt import BLOSUM62, translate

from ._aligner import Aligner


__all__ = ['validate']


def validate(refseq, seqs, score_matrix=None, mismatch=0, codon=True, revcomp=False, expected_identity=0., keep_insertions=True, quiet=False):
    msg = "cannot validate sequences that are not SeqRecord, Seq, or str objects"

    if isinstance(refseq, SeqRecord):
        r = str(refseq.seq)
    elif isinstance(refseq, Seq):
        r = str(refseq)
    elif isinstance(refseq, str):
        r = refseq
    else:
        raise ValueError(msg)

    qs = []
    for i, q in enumerate(seqs):
        if isinstance(q, SeqRecord):
            qs.append(str(q.seq))
        elif isinstance(q, Seq):
            qs.append(str(q))
        elif isinstance(q, str):
            qs.append(q)
        else:
            raise ValueError(msg)

    if score_matrix is None:
        score_matrix = BLOSUM62.load()

    aligner = Aligner(codon=codon)
    refs, queries, _, _, identities = aligner(
        r,
        qs,
        score_matrix,
        revcomp,
        expected_identity,
        keep_insertions,
        quiet
    )

    lengths = []
    scores = []
    for r, q, i in zip(refs, queries, identities):
        assert len(r) == len(q), 'sequences unaligned for some reason'
        lengths.append(len(r))
        if expected_identity > 0. and i < expected_identity:
            scores.append(None)
        else:
            scores.append(score_matrix(translate(r), translate(q), mismatch))

    return lengths, scores
