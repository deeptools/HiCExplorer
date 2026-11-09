"""Writes the matrix files of cpp/tests/cli/malformed_matrix_cli.sh.

Each file is a four bin matrix. good.cool and good.h5 are well formed; the
others carry one structural error that a reader must reject instead of
indexing past its arrays:

  bad_column.cool  a pixel whose bin2_id lies far outside the bin table
  bad_column.h5    a column index outside the matrix
  bad_indptr.h5    a last row offset past the stored entries

Usage: python make_malformed_fixtures.py <output directory>
(needs cooler, h5py and numpy; the files are committed, the test itself
needs no Python)
"""
import sys

import cooler
import h5py
import numpy as np
import pandas as pd

out = sys.argv[1]

bins = pd.DataFrame({"chrom": ["chr1"] * 4,
                     "start": [0, 10, 20, 30],
                     "end": [10, 20, 30, 40]})
pixels = pd.DataFrame({"bin1_id": [0, 0, 1, 2, 3],
                       "bin2_id": [0, 1, 1, 2, 3],
                       "count": [5, 2, 4, 3, 1]})
for name in ("good", "bad_column"):
    cooler.create_cooler(f"{out}/{name}.cool", bins, pixels,
                         assembly="test", h5opts={"compression": None})
with h5py.File(f"{out}/bad_column.cool", "r+") as f:
    f["pixels/bin2_id"][0] = 0x10000000


def write_h5(path, indptr, indices):
    with h5py.File(path, "w") as f:
        f.create_dataset("matrix/data", data=np.ones(len(indices)))
        f.create_dataset("matrix/indices", data=np.asarray(indices, dtype=np.int32))
        f.create_dataset("matrix/indptr", data=np.asarray(indptr, dtype=np.int64))
        f.create_dataset("matrix/shape", data=np.array([4, 4], dtype=np.int64))
        f.create_dataset("intervals/chr_list", data=np.array([b"chr1"] * 4))
        f.create_dataset("intervals/start_list", data=np.arange(0, 40, 10, dtype=np.int64))
        f.create_dataset("intervals/end_list", data=np.arange(10, 50, 10, dtype=np.int64))
        f.create_dataset("intervals/extra_list", data=np.ones(4))


write_h5(f"{out}/good.h5", [0, 2, 3, 4, 5], [0, 1, 1, 2, 3])
write_h5(f"{out}/bad_column.h5", [0, 2, 3, 4, 5], [0, 0x10000000, 1, 2, 3])
write_h5(f"{out}/bad_indptr.h5", [0, 2, 3, 4, 1005], [0, 1, 1, 2, 3])
