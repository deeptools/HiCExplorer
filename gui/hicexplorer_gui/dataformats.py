"""Detects the format of a data file a project loads.

The extension gives the first guess. Where the extension is ambiguous, the
content is checked the same way the C++ readers tell formats apart: an HDF5
signature can hold a HiCExplorer h5 matrix, a single-resolution cooler or a
multi-resolution mcool (told apart by their HDF5 group layout, as
cool_adapter.cpp and hic_adapter.cpp write them), a ``.hic`` file starts with
the ``HIC`` magic hicfilecpp's reader checks, a BAM file is a gzip (BGZF)
stream whose first block starts with ``BAM\\x01``, and a ``.pairs`` file
starts with the ``## pairs format v`` header pairs_file.cpp reads.

Returned formats use the same names the tool specifications use in a file
argument's ``formats`` list (``h5``, ``cool``, ``mcool``, ``hic``, ``bam``,
``sam``, ``pairs``, ``pairs.gz``), plus ``fastq``/``fastq.gz`` for the one
format no C++ tool reads directly. ``None`` means the format could not be
told, and the caller falls back to the unfiltered tool list.
"""

import gzip
import os

HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"
HIC_MAGIC = b"HIC"
GZIP_MAGIC = b"\x1f\x8b"
PAIRS_HEADER = b"## pairs format v"
BAM_MAGIC = b"BAM\x01"

FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
PAIRS_EXTENSIONS = (".pairs.gz", ".pairs")

# Every format a project can load, for file dialog filters.
LOADABLE_FORMATS = ("cool", "mcool", "h5", "hic", "bam", "sam", "pairs", "pairs.gz", "fastq", "fastq.gz")

FASTQ_MESSAGE = ("no available tool reads FASTQ directly; align it to BAM first, "
                 "then load the BAM file")


def _sniff_hdf5(path):
    """cool, mcool or h5 from the HDF5 group layout, or None when none matches."""
    try:
        import h5py
    except ImportError:
        return None
    try:
        with h5py.File(path, "r") as handle:
            keys = set(handle.keys())
    except OSError:
        return None
    if "resolutions" in keys:
        return "mcool"
    if "bins" in keys and "pixels" in keys:
        return "cool"
    if "intervals" in keys and "matrix" in keys:
        return "h5"
    return None


def _sniff_gzip_head(path, size=16):
    try:
        with gzip.open(path, "rb") as handle:
            return handle.read(size)
    except OSError:
        return b""


def _from_extension(lower):
    if lower.endswith(".mcool"):
        return "mcool"
    if lower.endswith(".cool"):
        return "cool"
    if lower.endswith(".h5"):
        return "h5"
    if lower.endswith(".hic"):
        return "hic"
    if lower.endswith(".bam"):
        return "bam"
    if lower.endswith(".sam"):
        return "sam"
    if lower.endswith(".pairs.gz"):
        return "pairs.gz"
    if lower.endswith(".pairs"):
        return "pairs"
    return None


def detect_format(path):
    """Best format guess for path, or None when nothing recognisable matches."""
    lower = path.lower()
    if any(lower.endswith(ext) for ext in FASTQ_EXTENSIONS):
        return "fastq.gz" if lower.endswith(".gz") else "fastq"
    try:
        with open(path, "rb") as handle:
            head = handle.read(4096)
    except OSError:
        return None
    if head[:8] == HDF5_SIGNATURE:
        return _sniff_hdf5(path) or _from_extension(lower) or "h5"
    if head[:3] == HIC_MAGIC:
        return "hic"
    if head[:2] == GZIP_MAGIC:
        block = _sniff_gzip_head(path)
        if block[:4] == BAM_MAGIC:
            return "bam"
        if block.startswith(PAIRS_HEADER):
            return "pairs.gz"
        return _from_extension(lower)
    if head.startswith(PAIRS_HEADER):
        return "pairs"
    return _from_extension(lower) or ("sam" if head[:1] == b"@" and lower.endswith(".sam") else None)


def is_matrix_format(fmt):
    return fmt in ("cool", "mcool", "h5", "hic")


def is_fastq_format(fmt):
    return fmt in ("fastq", "fastq.gz")
