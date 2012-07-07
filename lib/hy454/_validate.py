
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BioExt import (BLOSUM62,
    DNA80, translate)

from ._aligner import Aligner


__all__ = ['validate']


def validate(
    refseq,
    seqs,
    dna_score_matrix=None,
    protein_score_matrix=None,
    dna_mismatch=0,
    protein_mismatch=0,
    codon=True,
    revcomp=False,
    expected_identity=0.,
    keep_insertions=True,
    quiet=False):

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

    if dna_score_matrix is None:
        dna_score_matrix = DNA80

    if protein_score_matrix is None:
        score_matrix = BLOSUM62.load()

    if codon:
        score_matrix = protein_score_matrix
    else:
        score_matrix = dna_score_matrix

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
    dna_scores = []
    protein_scores = []
    for r, q, i in zip(refs, queries, identities):
        assert len(r) == len(q), 'sequences unaligned for some reason'
        lengths.append(len(r))
        if expected_identity > 0. and i < expected_identity:
            dna_scores.append(None)
            protein_scores.append(None)
        else:
            dna_scores.append(dna_score_matrix(r, q, dna_mismatch))
            protein_scores.append(
                protein_score_matrix(
                    translate(r),
                    translate(q),
                    protein_mismatch
                )
            )

    return lengths, dna_scores, protein_scores
