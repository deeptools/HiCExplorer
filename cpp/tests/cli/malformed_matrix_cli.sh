#!/usr/bin/env bash
# Matrix files whose sparse structure points outside the matrix are rejected
# with an error message and exit status 1. Before the structure check such a
# file made the tools index past their arrays.
#
#  1. A cool file with a pixel whose bin2_id lies far outside the bin table.
#  2. An h5 file with a column index outside the matrix.
#  3. An h5 file whose last row offset points past the stored entries.
#  4. The unmodified cool file and a well-formed h5 file of the same layout
#     still load.
#
# Usage: malformed_matrix_cli.sh <hicTransform binary> <hicInfo binary>
#                                <test_data directory> <python with h5py>

set -u
transform=$1
info=$2
data=$3
python=$4
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
failures=0

expect() {
    local description=$1
    shift
    if "$@"; then
        echo "ok: $description"
    else
        echo "FAILED: $description"
        failures=$((failures + 1))
    fi
}

# rejected <file> <label>: hicTransform exits 1 (not by a signal) and names the problem.
rejected() {
    local file=$1 label=$2
    "$transform" -m "$file" -o "$work/out_$label.cool" --method obs_exp \
        > /dev/null 2> "$work/err_$label.txt"
    local status=$?
    expect "$label: hicTransform exits 1, not by a signal (got $status)" [ "$status" -eq 1 ]
    expect "$label: the error names the invalid structure" \
        grep -q "not a valid" "$work/err_$label.txt"
}

cp "$data/small_test_matrix.cool" "$work/bad_column.cool"
chmod u+w "$work/bad_column.cool"

cat > "$work/make_files.py" <<'PYTHON'
import sys

import h5py
import numpy as np

work = sys.argv[1]

# The cool file: bin2_id is rewritten uncompressed with one pixel moved far
# outside the bin table.
with h5py.File(f"{work}/bad_column.cool", "r+") as f:
    bin2 = f["pixels/bin2_id"][()]
    bin2[0] = 0x10000000
    del f["pixels/bin2_id"]
    f.create_dataset("pixels/bin2_id", data=bin2)


def write_h5(path, indptr, indices):
    # A four bin HiCExplorer h5 matrix, uncompressed.
    with h5py.File(path, "w") as f:
        f.create_dataset("matrix/data", data=np.ones(len(indices)))
        f.create_dataset("matrix/indices", data=np.asarray(indices, dtype=np.int32))
        f.create_dataset("matrix/indptr", data=np.asarray(indptr, dtype=np.int64))
        f.create_dataset("matrix/shape", data=np.array([4, 4], dtype=np.int64))
        f.create_dataset("intervals/chr_list", data=np.array([b"chr1"] * 4))
        f.create_dataset("intervals/start_list", data=np.arange(0, 40, 10, dtype=np.int64))
        f.create_dataset("intervals/end_list", data=np.arange(10, 50, 10, dtype=np.int64))
        f.create_dataset("intervals/extra_list", data=np.ones(4))


write_h5(f"{work}/good.h5", [0, 2, 3, 4, 5], [0, 1, 1, 2, 3])
write_h5(f"{work}/bad_column.h5", [0, 2, 3, 4, 5], [0, 0x10000000, 1, 2, 3])
write_h5(f"{work}/bad_indptr.h5", [0, 2, 3, 4, 1005], [0, 1, 1, 2, 3])
PYTHON
"$python" "$work/make_files.py" "$work"
expect "the malformed files were written" [ $? -eq 0 ]

rejected "$work/bad_column.cool" cool_column
rejected "$work/bad_column.h5" h5_column
rejected "$work/bad_indptr.h5" h5_indptr

"$info" -m "$data/small_test_matrix.cool" > /dev/null 2>&1
expect "the unmodified cool file still loads" [ $? -eq 0 ]
"$transform" -m "$work/good.h5" -o "$work/good_out.h5" --method obs_exp > /dev/null 2>&1
expect "a well-formed h5 file of the same layout loads" [ $? -eq 0 ]

if [ "$failures" -ne 0 ]; then
    echo "$failures check(s) failed"
    exit 1
fi
echo "all checks passed"
