# -*- coding: utf-8 -*-
"""
Randomly split the nonredundant monthly sequence sets (2020-01 to 2024-02)
into 3 independent subsets, for robustness evaluation of the model
(comparison of results across independently sampled subsets of sequences).

Source data: data_preparation/monthly-cleanpd-2023/<YYYY-MM>.xlsx
Output: split_1/, split_2/, split_3/ (created in this folder), each containing
        one <YYYY-MM>.xlsx per month with a disjoint ~1/3 random subset of
        that month's rows.
"""
import os
import numpy as np
import pandas as pd

SEED = 42
N_SPLITS = 3

SOURCE_DIR = "/Volumes/Extreme SSD/data_preparation/monthly-cleanpd-2023"
DEST_DIR = os.path.dirname(os.path.abspath(__file__))

START_MONTH = "2020-01"
END_MONTH = "2024-02"


def month_range(start, end):
    start_year, start_mon = (int(x) for x in start.split("-"))
    end_year, end_mon = (int(x) for x in end.split("-"))
    months = []
    y, m = start_year, start_mon
    while (y, m) <= (end_year, end_mon):
        months.append(f"{y}-{m:02d}")
        m += 1
        if m > 12:
            m = 1
            y += 1
    return months


if __name__ == "__main__":
    rng = np.random.RandomState(SEED)

    split_dirs = [os.path.join(DEST_DIR, f"split_{i + 1}") for i in range(N_SPLITS)]
    for d in split_dirs:
        os.makedirs(d, exist_ok=True)

    for month in month_range(START_MONTH, END_MONTH):
        src_path = os.path.join(SOURCE_DIR, f"{month}.xlsx")
        if not os.path.exists(src_path):
            print(f"Skipping {month}: file not found at {src_path}")
            continue

        data = pd.read_excel(src_path)

        shuffled_idx = rng.permutation(len(data))
        chunks = np.array_split(shuffled_idx, N_SPLITS)

        for chunk_idx, split_dir in zip(chunks, split_dirs):
            subset = data.iloc[chunk_idx].sort_index()
            subset.to_excel(os.path.join(split_dir, f"{month}.xlsx"), index=False)

        print(f"{month}: {len(data)} rows -> {[len(c) for c in chunks]}")
