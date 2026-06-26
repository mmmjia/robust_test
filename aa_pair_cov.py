# -*- coding: utf-8 -*-
"""
For each month, find the top-m highest-covariance mutation pairs and track
how each pair's covariance value evolved across all preceding months.
Run independently on each of the 3 random splits (split_1/2/3) for the
robustness comparison requested in the review.

Input:  cov_results_split_<i>/covariance/<YYYY-MM>covariance.xlsx  (from cov_compu.py)
Output: cov_hilbert_split_<i>/<YYYY-MM>cov_hilbert.xlsx
"""
import os

import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SPLIT_IDS = [1, 2, 3]
M_TOP_PAIRS = 2000
START_MONTH_INDEX = 7  # skip the earliest, smallest months (too few mutations to be meaningful)

year = ['2020', '2021', '2022', '2023', '2024']
months = [str(i).zfill(2) for i in range(1, 12 + 1)]
month0 = [y + '-' + mo for y in year for mo in months]


def load_month_matrix(in_dir, month_str, cache):
    if month_str not in cache:
        df = pd.read_excel(os.path.join(in_dir, '{}covariance.xlsx'.format(month_str)), index_col=0, engine='calamine')
        labels = df.index.astype(str).tolist()
        label_to_idx = {lab: i for i, lab in enumerate(labels)}
        cache[month_str] = (labels, label_to_idx, df.values)
    return cache[month_str]


def run_split(split_id):
    in_dir = os.path.join(BASE_DIR, 'cov_results_split_{}'.format(split_id), 'covariance')
    out_dir = os.path.join(BASE_DIR, 'cov_hilbert_split_{}'.format(split_id))
    os.makedirs(out_dir, exist_ok=True)

    available = [m for m in month0 if os.path.exists(os.path.join(in_dir, '{}covariance.xlsx'.format(m)))]
    last_index = month0.index(available[-1])

    cache = {}
    for tt in range(START_MONTH_INDEX, last_index + 1):
        month_tt = month0[tt]
        print('split_{0} - {1}'.format(split_id, month_tt))

        siteaa, _, covariance = load_month_matrix(in_dir, month_tt, cache)
        covsele = np.triu(covariance)
        flattened_matrix = covsele.flatten()
        m = min(M_TOP_PAIRS, flattened_matrix.size)

        sorted_indices = np.argsort(flattened_matrix)[-m:]
        positions = np.unravel_index(sorted_indices, covsele.shape)

        variant_site1 = [siteaa[int(site)] for site in positions[0]]
        variant_site2 = [siteaa[int(site)] for site in positions[1]]

        cor_time = np.zeros((tt + 1, m))
        for i in range(m):
            cor_time[tt, i] = covariance[positions[0][i]][positions[1][i]]

        for mo in range(tt):
            siteaa_mo, idx_mo, covariance_mo = load_month_matrix(in_dir, month0[mo], cache)
            for i in range(m):
                indi = idx_mo.get(variant_site1[i])
                if indi is not None:
                    indj = idx_mo.get(variant_site2[i])
                    if indj is not None:
                        cor_time[mo, i] = covariance_mo[indi][indj]

        cor_time[cor_time < 0] = 0  # delete the negative covariance (the repulsion)
        df = pd.DataFrame(cor_time)
        df.columns = [str(variant_site1[i]) + '-' + str(variant_site2[i]) for i in range(m)]
        df.to_excel(os.path.join(out_dir, '{}cov_hilbert.xlsx'.format(month_tt)), index=False)


if __name__ == '__main__':
    for split_id in SPLIT_IDS:
        run_split(split_id)
