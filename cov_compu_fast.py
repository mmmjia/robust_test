# -*- coding: utf-8 -*-
"""
Vectorized (numpy/Counter) reimplementation of compu_mon_cov / covcaculate
from cov_compu.py. Mathematically equivalent to the original pandas-merge
based version, but avoids materializing the full O(n^2) mutation-pair
universe through pandas string joins, which is what made the original slow
for months with thousands of distinct mutations.

Kept as a separate module until validated against cov_compu.py's output.
"""

import math
import os
from collections import Counter

import numpy as np
import pandas as pd


def Unique_pair_func(mutList):
    return [mutList[i] + "|" + mutList[j] for i in range(len(mutList)) for j in range(len(mutList)) if i < j]


def _pos_of(mut):
    return int(''.join(ch for ch in mut if ch.isdigit()))


# Peak footprint of the fast path is dominated by a handful of n x n arrays
# (the float64 covariance matrix, the boolean same-position mask, and a
# transient copy while building the outer product) plus headroom for the
# DataFrame/xlsxwriter export. 24 bytes/cell is a conservative multiple of
# the 8-byte float64 matrix itself.
BYTES_PER_CELL_ESTIMATE = 24


def estimate_required_bytes(n):
    return BYTES_PER_CELL_ESTIMATE * n * n


def available_memory_bytes():
    try:
        import psutil
        return psutil.virtual_memory().available
    except ImportError:
        return None


def compu_mon_cov(dfInpMon, month_str, out_dir, save=True):
    mut_col = 'mutation' if 'mutation' in dfInpMon.columns else 'mutation info'
    Mutation0 = dfInpMon[mut_col].tolist()
    N = len(Mutation0)

    single_counter = Counter()
    pair_counter = Counter()
    for mutStr in Mutation0:
        mutList = mutStr.split(";")
        single_counter.update(mutList)
        pair_counter.update(Unique_pair_func(mutList))
    print('{0}: Total Seq = {1}'.format(month_str, N))

    # Same ordering rule as the original: by position ascending, count descending.
    MutIndex = sorted(single_counter.keys(), key=lambda m: (_pos_of(m), -single_counter[m]))
    idx = {m: i for i, m in enumerate(MutIndex)}
    pos_arr = np.array([_pos_of(m) for m in MutIndex])
    p = np.array([single_counter[m] for m in MutIndex], dtype=float) / N

    Cov_arr = -np.outer(p, p)  # baseline for non-co-occurring pairs: -P(x)P(y)
    same_pos = pos_arr[:, None] == pos_arr[None, :]
    Cov_arr[same_pos] = 0.0  # zeroes diagonal too, and same-residue (alt-AA) pairs

    for pair_str, count in pair_counter.items():
        a, b = pair_str.split('|')
        i, j = idx[a], idx[b]
        if pos_arr[i] == pos_arr[j]:
            continue
        cov_val = count / N - p[i] * p[j]
        Cov_arr[i, j] = cov_val
        Cov_arr[j, i] = cov_val

    if save:
        out_path = os.path.join(out_dir, '{}covariance.xlsx'.format(month_str))
        matrix_sele = pd.DataFrame(Cov_arr, columns=MutIndex, index=MutIndex)
        with pd.ExcelWriter(out_path, engine='xlsxwriter') as writer:
            matrix_sele.to_excel(writer, sheet_name='Sheet1')

    return Cov_arr, MutIndex


def covcaculate(covariance, sitew):  # calculate covariance at residue level
    siteonly = np.array([int(item.split('|')[0][1:-1]) for item in sitew])
    alpha = np.unique(siteonly)
    n = len(alpha)
    cov = np.zeros((n, n))

    site_indices = [np.where(siteonly == a)[0] for a in alpha]

    for s1 in range(n):
        u1 = site_indices[s1]
        for s2 in range(s1):
            u2 = site_indices[s2]
            cv2 = covariance[np.ix_(u1, u2)]
            pos = cv2[cv2 > 0]
            cov[s1, s2] = math.sqrt(np.sum(pos ** 2))
    cov = cov + cov.T
    return cov, alpha
