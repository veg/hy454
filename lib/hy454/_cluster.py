#!/usr/bin/env python3

import sys

from argparse import ArgumentParser, FileType
from collections import Counter
from json import dump as json_dump
from math import log
from re import compile as re_compile, I as re_I

import numpy as np

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import DBSCAN

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


DEBUG = True


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

def kmer_dist(k, m, a, b):
    ak = a.keys()
    bk = b.keys()
    s = np.sum((a[k] - b[k]) ** 2 for k in ak & bk)
    t = np.sum(a[k] ** 2 + b[k] ** 2 for k in ak ^ bk)
    return np.sqrt(s + t)

def embed(k, m, seqs):
    N = len(seqs)
    T = round(log(N) ** 2)
    step = round(N / T)
    # pull the reference sequences
    refs = list(set(str(s.seq) for i, s in enumerate(sorted(seqs, key=lambda s: -len(s))) if i % step == 0))
    T = len(refs)
    # get the counters for the seqs
    seq_counts = []
    for i, s in enumerate(seqs):
        seq_counts.append(kmer_counts(k, m, s))
        debug('computing %d-mer counts for %d seqs: %3d%%\r' % (k, N, 100 * (i + 1) // N), end='')
    debug()
    # get the counters for the refs
    ref_counts = []
    for i, r in enumerate(refs):
        ref_counts.append(kmer_counts(k, m, r))
        debug('computing %d-mer counts for %d refs: %3d%%\r' % (k, T, 100 * (i + 1) // T), end='')
    debug()
    # populate the reference distance matrix
    Q = np.zeros((T, T), dtype=float)
    for i, s in enumerate(ref_counts):
        start = i + 1
        for j, r in enumerate(ref_counts[start:], start=start):
            d = kmer_dist(k, m, s, r)
            Q[i, j] = d
            Q[j, i] = d
            debug('embedding %d references in %d dims: %3d%%\r' % (T, T, (100 * (i * T + j + 1) // (T * T), end='')
    debug()
    # populate the embedding matrix
    R = np.zeros((N, T), dtype=float)
    ncomp = N * T
    for i, s in enumerate(seq_counts):
        for j, r in enumerate(ref_counts):
            R[i, j] = kmer_dist(k, m, s, r)
            debug('embedding %d sequences in %d dims: %3d%%\r' % (N, T, 100 * (i * T + j + 1) // ncomp), end='')
    debug()
    # normalize to unit lengths
    d = np.seterr(all='ignore')
    R = np.nan_to_num(R / (np.sum(R, axis=1)[:, np.newaxis]))
    np.seterr(**d)
    return Q, R

def cluster(seqs, k, m):
    Q, R = embed(k, m, seqs)
    # debug('computing pairwise distances')
    D = pairwise_distances(R, metric='cosine', n_jobs=-1)
    debug('clustering using dbscan')
    db = DBSCAN(eps=0.05, min_samples=1, metric='precomputed').fit(D)
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
