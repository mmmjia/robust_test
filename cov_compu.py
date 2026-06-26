# -*- coding: utf-8 -*-
"""
Per-month mutation-pair and residue-level covariance, computed independently
on each of the 3 random splits (split_1, split_2, split_3) created by
split_monthly_data.py, for the robustness comparison requested in the review.

For each split_<i>/<YYYY-MM>.xlsx, this writes:
  cov_results_split_<i>/covariance/<YYYY-MM>covariance.xlsx          (mutation-pair level)
  cov_results_split_<i>/covariance_residue_level/<YYYY-MM>.xlsx      (residue level)
"""

import math
import os

import numpy as np
import pandas as pd

import ImportantFunc as Imp
import cov_compu_fast as fast

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SPLIT_IDS = [1, 2, 3]
MEM_SAFETY_FRACTION = 0.5  # only take the fast path if it needs less than this share of currently free RAM


def Unique_pair_func(mutList):
    return [mutList[i] + "|" + mutList[j] for i in range(len(mutList)) for j in range(len(mutList)) if i < j]


def compu_mon_cov(dfInpMon, month_str, out_dir, save=True):
    pd.set_option('mode.chained_assignment', None)  # deactivate the copywarning from pandas

    mut_col = 'mutation' if 'mutation' in dfInpMon.columns else 'mutation info'
    Mutation0 = dfInpMon[mut_col].tolist()
    MutPair, MutSingle = [], []  # Mutation0: list of per-sequence mutation strings
    for mutList in Mutation0:
        mutList = mutList.split(";")
        MutSingle.extend(mutList)
        MutPair.extend(Unique_pair_func(mutList))
    print('{0}: Total Seq = {1}'.format(month_str, len(Mutation0)))

    M = Imp.count_dups(MutPair); dfCountPair = pd.DataFrame({'Mut-Pair': M[0], 'CountPair': M[1]})
    n = Imp.count_dups(MutSingle); dfSingleX = pd.DataFrame({'x': n[0], 'CountX': n[1]})
    dfSingleY = pd.DataFrame({'y': n[0], 'CountY': n[1]})
    MutEle = [ele for Mut in Mutation0 for ele in Mut.split(";")]
    M = Imp.count_dups(sorted(MutEle)); df = pd.DataFrame({'Mutation': M[0], 'Count': M[1]})
    df['Pos'] = df['Mutation'].astype(str).str.extractall(r'(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    df.sort_values(['Pos', 'Count'], inplace=True, ascending=[True, False])
    df = df.reset_index(drop=True)
    MutIndex = df['Mutation'].tolist()  # MutIndex is the list of sorted distinct mutations

    dfMutPair = pd.DataFrame({'Mut-Pair': Unique_pair_func(MutIndex)})
    dfMutPair[['x', 'y']] = dfMutPair['Mut-Pair'].str.split('|', n=1, expand=True)
    dfMutPair = pd.merge(dfMutPair, dfCountPair, how='left', on='Mut-Pair')
    dfMutPair = pd.merge(dfMutPair, dfSingleX, how='left', on='x')
    dfMutPair = pd.merge(dfMutPair, dfSingleY, how='left', on='y').fillna(0)
    dfMutPair['Covariant'] = (dfMutPair['CountPair'] / len(Mutation0)) - (dfMutPair['CountX'] / len(Mutation0)) * (dfMutPair['CountY'] / len(Mutation0))
    dfMutPair['PosX'] = dfMutPair['x'].astype(str).str.extractall(r'(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    dfMutPair['PosY'] = dfMutPair['y'].astype(str).str.extractall(r'(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    dfMutPair['Marker'] = dfMutPair.apply(lambda x: 1 if x['PosX'] == x['PosY'] else 0, axis=1)
    dfMutPair.loc[dfMutPair['Marker'] == 1, 'Covariant'] = 0

    dfMutPair1 = dfMutPair[['x', 'y', 'Covariant']].copy(deep=True)
    dfMutPair1['Mut-Pair'] = dfMutPair['x'] + "|" + dfMutPair['y']; dfMutPair1.drop(['x', 'y'], inplace=True, axis=1)
    dfMutPair2 = dfMutPair[['y', 'x', 'Covariant']].copy(deep=True)
    dfMutPair2['Mut-Pair'] = dfMutPair['y'] + "|" + dfMutPair['x']; dfMutPair2.drop(['x', 'y'], inplace=True, axis=1)
    dfMutPair = pd.concat([dfMutPair1, dfMutPair2], axis=0).reset_index(drop=True); dfMutPair = dfMutPair[['Mut-Pair', 'Covariant']]
    dfMutFrame = pd.DataFrame({'Mut-Pair': [MutIndex[i] + "|" + MutIndex[j] for i in range(len(MutIndex)) for j in range(len(MutIndex))]})
    dfMutFrame = pd.merge(dfMutFrame, dfMutPair, how='left', on='Mut-Pair').fillna(0)
    Cov_arr = np.reshape(dfMutFrame['Covariant'].to_numpy(), (len(MutIndex), len(MutIndex)))

    if save:
        out_path = os.path.join(out_dir, '{}covariance.xlsx'.format(month_str))
        matrix_sele = pd.DataFrame(Cov_arr, columns=MutIndex, index=MutIndex)
        with pd.ExcelWriter(out_path, engine='xlsxwriter') as writer:
            matrix_sele.to_excel(writer, sheet_name='Sheet1')

    return Cov_arr, MutIndex


def covcaculate(covariance, sitew):  # calculate covariance at residue level
    siteonly = []
    for item in sitew:
        item = item.split('|')[0]
        site0 = item[1:-1]
        siteonly.append(int(site0))

    alpha = np.unique(siteonly)
    n = np.size(alpha)
    cov = np.zeros((n, n))
    nn = len(siteonly)

    class gresks:
        def __init__(self, x, y):
            self.x = x
            self.y = y

    sites = []
    for i in range(n):
        ii = 0
        indi = []
        while ii < nn and int(siteonly[ii]) < int(alpha[i]) + 1:
            if siteonly[ii] == alpha[i]:
                indi.append(ii)
            ii = ii + 1
        sites.append(gresks(alpha[i], indi))
    for s1 in range(n):
        posall = sites[s1]
        u1 = posall.y
        for s2 in range(s1):
            posall2 = sites[s2]
            u2 = posall2.y
            cv2 = covariance[u1, :][:, u2]
            sum_squares = 0
            for row in cv2:
                for element in row:
                    if element > 0:
                        sum_squares += element ** 2
            cov[s1, s2] = math.sqrt(sum_squares)
    cov = cov + cov.T
    return cov, alpha


def count_unique_mutations(dfInpMon):
    mut_col = 'mutation' if 'mutation' in dfInpMon.columns else 'mutation info'
    muts = set()
    for m in dfInpMon[mut_col]:
        muts.update(m.split(';'))
    return len(muts)


def compu_mon_cov_checked(dfInpMon, month_str, out_dir, save=True):
    """Use the fast path only if it fits comfortably in currently available RAM;
    otherwise (or if it still hits MemoryError) fall back to the slow, validated
    original implementation in this module."""
    n = count_unique_mutations(dfInpMon)
    need = fast.estimate_required_bytes(n)
    avail = fast.available_memory_bytes()

    if avail is not None and need < avail * MEM_SAFETY_FRACTION:
        try:
            Cov_arr, MutIndex = fast.compu_mon_cov(dfInpMon, month_str, out_dir, save=save)
            return Cov_arr, MutIndex, 'fast'
        except MemoryError:
            print('{0}: fast path hit MemoryError (n={1}), falling back to original implementation'.format(month_str, n))
    else:
        avail_mb = 'unknown' if avail is None else '{:.0f}MB'.format(avail / 1e6)
        print('{0}: n={1}, need~{2:.0f}MB > safety threshold (avail={3}), using original implementation'.format(
            month_str, n, need / 1e6, avail_mb))

    Cov_arr, MutIndex = compu_mon_cov(dfInpMon, month_str, out_dir, save=save)
    return Cov_arr, MutIndex, 'slow'


def covcaculate_checked(covariance, sitew):
    n = len(sitew)
    need = fast.estimate_required_bytes(n)
    avail = fast.available_memory_bytes()

    if avail is not None and need < avail * MEM_SAFETY_FRACTION:
        try:
            return fast.covcaculate(covariance, sitew)
        except MemoryError:
            print('covcaculate: fast path hit MemoryError (n={0}), falling back to original implementation'.format(n))

    return covcaculate(covariance, sitew)


def run_split(split_id):
    split_dir = os.path.join(BASE_DIR, 'split_{}'.format(split_id))
    out_root = os.path.join(BASE_DIR, 'cov_results_split_{}'.format(split_id))
    out_pair_dir = os.path.join(out_root, 'covariance')
    out_residue_dir = os.path.join(out_root, 'covariance_residue_level')
    os.makedirs(out_pair_dir, exist_ok=True)
    os.makedirs(out_residue_dir, exist_ok=True)

    months = sorted(f[:-5] for f in os.listdir(split_dir) if f.endswith('.xlsx'))

    for month_str in months:
        print('split_{0} - {1}'.format(split_id, month_str))
        dfInp = pd.read_excel(os.path.join(split_dir, '{}.xlsx'.format(month_str)))
        Cov_arr, MutIndex, path_used = compu_mon_cov_checked(dfInp, month_str, out_pair_dir, save=True)
        print('  -> compu_mon_cov used {} path'.format(path_used))
        cov, alpha = covcaculate_checked(Cov_arr, MutIndex)
        matrix_sele = pd.DataFrame(cov, columns=alpha, index=alpha)
        matrix_sele.to_excel(os.path.join(out_residue_dir, '{}.xlsx'.format(month_str)))


if __name__ == '__main__':
    for split_id in SPLIT_IDS:
        run_split(split_id)
