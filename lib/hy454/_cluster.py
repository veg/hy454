#!/usr/bin/env python3

import sys

from argparse import ArgumentParser, FileType
from collections import Counter
from json import dump as json_dump
from math import log
from re import compile as re_compile, I as re_I

import numpy as np

from scipy.spatial.distance import cosine

from sklearn.cluster import DBSCAN
from sklearn.externals.joblib import Parallel, delayed
from sklearn.metrics.pairwise import pairwise_distances

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


DEBUG = True
N_JOBS = -1
PRE_DISPATCH = '5 * n_jobs'
VERBOSE = 3

def debug(string='', *args, **kwargs):
    if DEBUG:
        print(string, file=sys.stderr, *args, **kwargs)

hetero = re_compile(r'([a-z])\1+', re_I)
def kmer_counts(k, m, s):
    if isinstance(s, SeqRecord):
        s = str(s.seq)
    elif isinstance(s, Seq):
        s = str(s)
    elif not isinstance(s, str):
        raise ValueError('inputs must be SeqRecord, Seq, or str')
    s = s.replace('-', '')
    c = Counter()
    for lwr in range(0, len(s)):
        upr = min(lwr + k, len(s))
        p = hetero.sub(r'\1', s[lwr:upr])
        if len(p) >= m:
            c.update(p)
    return c

def kmer_dist(a, b):
    ak = a.keys()
    bk = b.keys()
    s = np.sum((a[k] - b[k]) ** 2 for k in ak & bk)
    t = np.sum(a[k] ** 2 + b[k] ** 2 for k in ak ^ bk)
    return np.sqrt(s + t)

def embed(k, m, seqs, dispatcher=None):
    if dispatcher is None:
        dispatcher = Parallel(n_jobs=N_JOBS, verbose=VERBOSE, pre_dispatch=PRE_DISPATCH)
    N = len(seqs)
    T = round(log(N) ** 2)
    step = round(N / T)
    # pull the reference sequences
    refs = list(set(str(s.seq) for i, s in enumerate(sorted(seqs, key=lambda s: -len(s))) if i % step == 0))
    T = len(refs)
    # get the counters for the seqs
    seq_counts = dispatcher(delayed(kmer_counts)(k, m, s) for s in seqs)
    assert len(seq_counts) == N
    # get the counters for the refs
    ref_counts = dispatcher(delayed(kmer_counts)(k, m, r) for r in refs)
    assert len(ref_counts) == T
    # populate the reference distance matrix
    Q = np.zeros((T, T), dtype=float)
    debug('computing distance matrix for %d references' % T)
    qs = dispatcher(delayed(kmer_dist)(s, r) for i, s in enumerate(ref_counts) for r in ref_counts[i + 1:])
    assert len(qs) == (T * (T - 1) // 2)
    idx = 0
    for i in range(T):
        for j in range(i + 1, T):
            Q[i, j] = qs[idx]
            Q[j, i] = qs[idx]
            idx += 1
    # populate the embedding matrix
    debug('embedding %d sequences in %d dimensions' % (N, T))
    rs = dispatcher(delayed(kmer_dist)(s, r) for s in seq_counts for r in ref_counts)
    assert len(rs) == (N * T)
    R = np.array(rs, dtype=float).reshape((N, T))
    # normalize to unit lengths
    d = np.seterr(all='ignore')
    R = np.nan_to_num(R / (np.sum(R, axis=1)[:, np.newaxis]))
    np.seterr(**d)
    return Q, R

def cluster(seqs, k, m):
    dispatcher = Parallel(n_jobs=N_JOBS, verbose=VERBOSE, pre_dispatch=PRE_DISPATCH)
    Q, R = embed(k, m, seqs, dispatcher)
    N, T = R.shape
    debug('computing pairwise distances')
    ds = dispatcher(delayed(cosine)(R[i, :], R[j, :]) for i in range(N) for j in range(i + 1, N))
    D = np.zeros((N, N), dtype=float)
    idx = 0
    for i in range(N):
        for j in range(i + 1, N):
            D[i, j] = ds[idx]
            D[j, i] = ds[idx]
            idx += 1
    # cluster
    debug('clustering using dbscan')
    db = DBSCAN(eps=0.01, min_samples=10, metric='precomputed').fit(D)
    # core = db.core_sample_indices_
    labels = db.labels_
    labelset = set(labels)
    n_clusters_ = len(labelset) - (1 if -1 in labelset else 0)
    debug('dbscan: found %d clusters' % n_clusters_)
    m = {}
    for l in labelset:
        m[l] = []
        for i, s in enumerate(seqs):
            if labels[i] == l:
                m[l].append(str(s.seq))
    return m

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = ArgumentParser(description='cluster sequences using DBSCAN')
    parser.add_argument('input', nargs='?', type=FileType('r'), default=sys.stdin)
    parser.add_argument('output', nargs='?', type=FileType('w'), default=sys.stdout)
    parser.add_argument('-k', type=int, default=16, help='K-mer size for hashing')
    parser.add_argument('-m', type=int, default=8, help='min size of K-mer after homopolymer reduction')

    ns = parser.parse_args(args)

    seqs = [r for r in SeqIO.parse(ns.input, 'fasta')]

    m = cluster(seqs, ns.k, ns.m)

    json_dump(m, sys.stdout, indent=1)

    return 0


if __name__ == '__main__':
    sys.exit(main())
