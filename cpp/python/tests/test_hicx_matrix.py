"""hicx_matrix against the reference readers, E2 (cpp/PLAN.md tier 10, item 10.2).

cool and mcool against cooler's Cooler.matrix(balance=...).fetch, h5 against
hicmatrix's hiCMatrix sliced with getRegionBinRange, .hic against hicstraw
1.3.1's getRecordsAsMatrix. See support.py for the environment variables.
"""

import os

import numpy as np
import pytest

from support import TEST_DATA, assert_e2, data, import_hicx_matrix, run_oracle

hicx_matrix = import_hicx_matrix()

GM = "hicTADClassifier/gm12878_chr1.cool"
KR = "hicCorrectMatrix/gm12878_KR.cool"
LOOPS = "hicDetectLoops/GSE63525_GM12878_insitu_primary_2_5mb.cool"
MCOOL = "hicBuildMatrix/multi_small_test_matrix.mcool"
GSM_H5 = "hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1.h5"
GSM_COOL = "hicDifferentialTAD/GSM2644945_Untreated-R1.100000_chr1_chr2.cool"
LI_H5 = "Li_et_al_2015.h5"
SMALL_HIC = "hicHyperoptDetectLoopsHiCCUPS/SRR1791297_30.hic"
LARGE_HIC = os.environ.get("HICX_LARGE_HIC", "")


def need(path):
    if not os.path.isfile(path.split("::")[0]):
        pytest.skip(f"{path} is not available")


# --------------------------------------------------------------------- cool, mcool

# id, file, URI suffix given to open(), resolution, region1, region2, normalization
COOL_CASES = [
    ("gm12878_diag_2mb", GM, "", None, "1:10,000,000-12,000,000", None, "none"),
    ("gm12878_offdiag", GM, "", None, "1:10,000,000-12,000,000", "1:20,000,000-21,500,000", "none"),
    ("gm12878_offdiag_lower", GM, "", 10000, "1:20,000,000-21,500,000", "1:10,000,000-12,000,000", "none"),
    ("gm12878_start", GM, "", None, "1:0-1,500,000", None, "none"),
    ("gm12878_end", GM, "", None, "1:248,000,000-249,250,621", None, "none"),
    ("gm12878_unaligned", GM, "", None, "1:10,005,000-10,995,001", "1:10,000,001-10,600,000", "none"),
    ("kr_raw_chr1", KR, "", None, "1", None, "none"),
    ("kr_balanced_chr1", KR, "", None, "1", None, "balanced"),
    ("kr_balanced_inter", KR, "", None, "1:0-50,000,000", "2", "balanced"),
    ("loops_KR", LOOPS, "", None, "3", None, "KR"),
    ("loops_VC_inter", LOOPS, "", None, "1:100,000,000-150,000,000", "2", "VC"),
    ("loops_VC_SQRT", LOOPS, "", None, "2:0-60,000,000", None, "VC_SQRT"),
    ("mcool_10kb_diag", MCOOL, "", 10000, "chr2L:1,000,000-3,000,000", None, "none"),
    ("mcool_10kb_chrom", MCOOL, "", 10000, "chr4", None, "none"),
    ("mcool_20kb_inter", MCOOL, "", 20000, "chr2L:0-1,000,000", "chr3R:5,000,000-6,000,000", "none"),
    ("mcool_20kb_chrom_end", MCOOL, "", 20000, "chrX:21,000,000-22,420,000", None, "none"),
    ("mcool_uri_5kb", MCOOL, "::/resolutions/5000", None, "chr3L:2,000,000-3,000,000",
     "chr3L:2,500,000-4,000,000", "none"),
    # 'format' stored as a fixed length string: cooler.fileops.is_cooler says
    # no, cooler.Cooler opens it, and so must hicx_matrix.open.
    ("gsm_fixed_format_chr1", GSM_COOL, "", None, "chr1:110,000,000-120,000,000", None, "none"),
    ("gsm_fixed_format_inter", GSM_COOL, "", None, "chr1:0-40,000,000", "chr2", "none"),
    ("gsm_fixed_format_chrom", GSM_COOL, "", 100000, "chr2", None, "none"),
]
COOL_BY_ID = {case[0]: case for case in COOL_CASES}


def cool_uri(case):
    _, relative, suffix, resolution, *_ = case
    if suffix:
        return data(relative) + suffix
    if relative == MCOOL:
        return f"{data(relative)}::/resolutions/{resolution}"
    return data(relative)


def cool_balance(normalization):
    return {"none": False, "balanced": True}.get(normalization, normalization)


@pytest.fixture(scope="module")
def cooler_ref(tmp_path_factory):
    queries = []
    for case in COOL_CASES:
        if os.path.isfile(data(case[1])):
            queries.append({"id": case[0], "what": "matrix", "uri": cool_uri(case),
                            "region1": case[4], "region2": case[5],
                            "balance": cool_balance(case[6])})
    for name, relative in (("gm12878", GM), ("mcool_5kb", MCOOL + "::/resolutions/5000"),
                           ("gsm_fixed_format", GSM_COOL)):
        if os.path.isfile(data(relative.split("::")[0])):
            queries.append({"id": name, "what": "chromosomes", "uri": data(relative)})
    return run_oracle("cooler", queries, tmp_path_factory.mktemp("cooler"))


@pytest.mark.parametrize("case", COOL_CASES, ids=[case[0] for case in COOL_CASES])
def test_cool_region(case, cooler_ref):
    case_id, relative, suffix, resolution, region1, region2, normalization = case
    need(data(relative))
    matrix = hicx_matrix.open(data(relative) + suffix)
    actual = matrix.fetch(region1, region2, resolution=resolution, normalization=normalization)
    assert_e2(actual, cooler_ref[case_id]["matrix"], "cooler")


@pytest.mark.parametrize("transform", ["log1p", "log"])
@pytest.mark.parametrize("case_id", ["gm12878_offdiag", "kr_balanced_chr1"])
def test_cool_transform(case_id, transform, cooler_ref):
    _, relative, suffix, resolution, region1, region2, normalization = COOL_BY_ID[case_id]
    need(data(relative))
    matrix = hicx_matrix.open(data(relative) + suffix)
    actual = matrix.fetch(region1, region2, resolution=resolution, normalization=normalization,
                          transform=transform)
    with np.errstate(divide="ignore", invalid="ignore"):
        expected = getattr(np, transform)(cooler_ref[case_id]["matrix"])
    assert_e2(actual, expected, "cooler")


def test_cool_metadata(cooler_ref):
    need(data(GM))
    gm = hicx_matrix.open(data(GM))
    assert gm.format == "cool"
    assert gm.resolutions() == [10000]
    assert gm.normalizations() == ["none"]
    reference = cooler_ref["gm12878"]
    assert gm.chromosomes() == list(zip(reference["names"].tolist(), reference["lengths"].tolist()))
    assert hicx_matrix.open(data(KR)).normalizations() == ["none", "balanced"]
    assert hicx_matrix.open(data(LOOPS)).normalizations() == ["none", "KR", "VC", "VC_SQRT"]


def test_cool_with_fixed_length_format_attribute(cooler_ref):
    need(data(GSM_COOL))
    matrix = hicx_matrix.open(data(GSM_COOL))
    assert matrix.format == "cool"
    assert matrix.resolutions() == [100000]
    reference = cooler_ref["gsm_fixed_format"]
    assert matrix.chromosomes() == list(zip(reference["names"].tolist(), reference["lengths"].tolist()))


def test_mcool_metadata(cooler_ref):
    need(data(MCOOL))
    mcool = hicx_matrix.open(data(MCOOL))
    assert mcool.format == "mcool"
    assert mcool.resolutions() == [5000, 10000, 20000]
    assert mcool.normalizations() == ["none"]
    reference = cooler_ref["mcool_5kb"]
    assert mcool.chromosomes() == list(zip(reference["names"].tolist(), reference["lengths"].tolist()))
    single = hicx_matrix.open(data(MCOOL) + "::/resolutions/20000")
    assert single.format == "cool"
    assert single.resolutions() == [20000]


# --------------------------------------------------------------------- h5

# id, file, region1 (chrom, start, end), region2; start and end None mean the
# whole chromosome.
H5_CASES = [
    ("gsm_diag_20mb", GSM_H5, ("chr1", 10_000_000, 30_000_000), None),
    ("gsm_offdiag", GSM_H5, ("chr1", 10_000_000, 20_000_000), ("chr1", 50_000_000, 60_000_000)),
    ("gsm_offdiag_lower", GSM_H5, ("chr1", 50_000_000, 60_000_000), ("chr1", 10_000_000, 20_000_000)),
    ("gsm_overlapping", GSM_H5, ("chr1", 15_050_000, 25_000_000), ("chr1", 20_000_000, 40_000_001)),
    ("gsm_start", GSM_H5, ("chr1", 0, 5_000_000), None),
    ("gsm_end", GSM_H5, ("chr1", 190_000_000, 197_195_432), None),
    ("gsm_chromosome", GSM_H5, ("chr1", None, None), None),
    ("li_fragments", LI_H5, ("X", 1_000_000, 1_500_000), ("X", 1_200_000, 2_000_123)),
    ("li_fragments_lower", LI_H5, ("X", 5_000_000, 5_400_000), ("X", 4_000_000, 4_300_000)),
    ("li_end", LI_H5, ("X", 21_500_000, None), ("X", 20_000_000, None)),
]


def region_text(region, lengths):
    chrom, start, end = region
    if start is None and end is None:
        return chrom
    end = lengths[chrom] if end is None else end
    return f"{chrom}:{start:,}-{end:,}"


@pytest.fixture(scope="module")
def hicmatrix_ref(tmp_path_factory):
    queries = []
    for case_id, relative, region1, region2 in H5_CASES:
        if os.path.isfile(data(relative)):
            queries.append({"id": case_id, "what": "matrix", "path": data(relative),
                            "region1": list(region1), "region2": list(region2 or region1)})
    for name, relative in (("gsm", GSM_H5), ("li", LI_H5)):
        if os.path.isfile(data(relative)):
            queries.append({"id": name, "what": "chromosomes", "path": data(relative)})
    return run_oracle("hicmatrix", queries, tmp_path_factory.mktemp("hicmatrix"))


@pytest.mark.parametrize("case", H5_CASES, ids=[case[0] for case in H5_CASES])
def test_h5_region(case, hicmatrix_ref):
    case_id, relative, region1, region2 = case
    need(data(relative))
    matrix = hicx_matrix.open(data(relative))
    lengths = dict(matrix.chromosomes())
    actual = matrix.fetch(region_text(region1, lengths),
                          None if region2 is None else region_text(region2, lengths))
    assert_e2(actual, hicmatrix_ref[case_id]["matrix"], "hicmatrix")


def test_h5_metadata(hicmatrix_ref):
    need(data(GSM_H5))
    gsm = hicx_matrix.open(data(GSM_H5))
    assert gsm.format == "h5"
    assert gsm.resolutions() == [100000]
    assert gsm.normalizations() == ["none"]
    for name, relative in (("gsm", GSM_H5), ("li", LI_H5)):
        reference = hicmatrix_ref[name]
        assert hicx_matrix.open(data(relative)).chromosomes() == \
            list(zip(reference["names"].tolist(), reference["lengths"].tolist()))
    # Restriction fragment bins have no single bin size.
    assert hicx_matrix.open(data(LI_H5)).resolutions() == []


# --------------------------------------------------------------------- .hic

SMALL_HIC_REGIONS = [
    ("intra", ("NC_001136.10", 100_000, 600_000), None),
    ("inter_swapped", ("NC_001136.10", 0, 500_000), ("NC_001133.9", None, None)),
    ("inter", ("NC_001133.9", 50_000, 200_000), ("NC_001136.10", 1_000_000, 1_531_933)),
    ("unaligned", ("NC_001136.10", 102_345, 297_001), ("NC_001136.10", 250_001, 400_000)),
    ("chromosome", ("NC_001133.9", None, None), None),
]
SMALL_HIC_SETTINGS = [(resolution, normalization, transform)
                      for resolution in (5000, 10000)
                      for normalization, transform in (("none", "none"), ("KR", "none"), ("VC", "none"),
                                                       ("none", "obs_exp"), ("KR", "obs_exp"))]

LARGE_HIC_CASES = [
    ("diag_2mb", 10000, "none", "none", ("1", 50_000_000, 52_000_000), None),
    ("diag_2mb", 10000, "KR", "none", ("1", 50_000_000, 52_000_000), None),
    ("diag_2mb", 25000, "none", "none", ("1", 50_000_000, 52_000_000), None),
    ("diag_2mb", 25000, "KR", "none", ("1", 50_000_000, 52_000_000), None),
    ("diag_2mb", 25000, "KR", "obs_exp", ("1", 50_000_000, 52_000_000), None),
    ("inter_swapped", 25000, "none", "none", ("2", 30_000_000, 31_000_000), ("1", 50_000_000, 51_000_000)),
    ("inter_swapped", 25000, "KR", "none", ("2", 30_000_000, 31_000_000), ("1", 50_000_000, 51_000_000)),
    ("chr1_start_without_records", 10000, "KR", "none", ("1", 0, 2_000_000), None),
]


def hic_case(prefix, path, name, resolution, normalization, transform, region1, region2):
    case_id = f"{prefix}_{name}_{resolution}_{normalization}_{transform}"
    return (case_id, path, resolution, normalization, transform, region1, region2)


SMALL_HIC_CASES = [hic_case("small", data(SMALL_HIC), name, *setting, region1, region2)
                   for setting in SMALL_HIC_SETTINGS for name, region1, region2 in SMALL_HIC_REGIONS]
LARGE_CASES = [hic_case("large", LARGE_HIC, *case) for case in LARGE_HIC_CASES]


def hicstraw_queries(cases):
    return [{"id": case_id, "what": "matrix", "path": path, "resolution": resolution,
             "norm": "NONE" if normalization == "none" else normalization,
             "matrix_type": "oe" if transform == "obs_exp" else "observed",
             "region1": list(region1), "region2": list(region2 or region1)}
            for case_id, path, resolution, normalization, transform, region1, region2 in cases]


@pytest.fixture(scope="module")
def small_hic_ref(tmp_path_factory):
    need(data(SMALL_HIC))
    queries = hicstraw_queries(SMALL_HIC_CASES)
    queries.append({"id": "small", "what": "chromosomes", "path": data(SMALL_HIC)})
    return run_oracle("hicstraw", queries, tmp_path_factory.mktemp("hicstraw_small"))


@pytest.fixture(scope="module")
def large_hic_ref(tmp_path_factory):
    if not LARGE_HIC or not os.path.isfile(LARGE_HIC):
        pytest.skip("HICX_LARGE_HIC does not name an existing .hic file")
    return run_oracle("hicstraw", hicstraw_queries(LARGE_CASES), tmp_path_factory.mktemp("hicstraw_large"))


def check_hic_case(case, reference):
    case_id, path, resolution, normalization, transform, region1, region2 = case
    matrix = hicx_matrix.open(path)
    lengths = dict(matrix.chromosomes())
    actual = matrix.fetch(region_text(region1, lengths),
                          None if region2 is None else region_text(region2, lengths),
                          resolution=resolution, normalization=normalization, transform=transform)
    assert_e2(actual, reference[case_id]["matrix"], "hicstraw")


@pytest.mark.parametrize("case", SMALL_HIC_CASES, ids=[case[0] for case in SMALL_HIC_CASES])
def test_hic_region(case, small_hic_ref):
    check_hic_case(case, small_hic_ref)


@pytest.mark.parametrize("case", LARGE_CASES, ids=[case[0] for case in LARGE_CASES])
def test_large_hic_region(case, large_hic_ref):
    check_hic_case(case, large_hic_ref)


def test_hic_metadata(small_hic_ref):
    hic = hicx_matrix.open(data(SMALL_HIC))
    reference = small_hic_ref["small"]
    assert hic.format == "hic"
    assert hic.chromosomes() == list(zip(reference["names"].tolist(), reference["lengths"].tolist()))
    assert hic.resolutions() == sorted(reference["resolutions"].tolist())
    norms = hic.normalizations()
    assert norms[0] == "none" and "KR" in norms and "VC" in norms and "NONE" not in norms


# --------------------------------------------------------------------- errors

def test_errors_cool():
    need(data(GM))
    gm = hicx_matrix.open(data(GM))
    with pytest.raises(KeyError, match="unknown chromosome 'chr99'"):
        gm.fetch("chr99:0-100")
    with pytest.raises(ValueError, match="malformed region|invalid region"):
        gm.fetch("1:abc-def")
    with pytest.raises(ValueError, match="region '1:5,000,000-1,000,000'"):
        gm.fetch("1:5,000,000-1,000,000")
    with pytest.raises(ValueError, match="invalid region"):
        gm.fetch("1:249,000,000-250,000,000")
    with pytest.raises(ValueError, match="resolution 5000"):
        gm.fetch("1:0-1,000,000", resolution=5000)
    with pytest.raises(ValueError, match="unknown normalization 'balanced'"):
        gm.fetch("1:0-1,000,000", normalization="balanced")
    with pytest.raises(ValueError, match="obs_exp"):
        gm.fetch("1:0-1,000,000", transform="obs_exp")
    with pytest.raises(ValueError, match="unknown transform"):
        gm.fetch("1:0-1,000,000", transform="sqrt")
    mcool = hicx_matrix.open(data(MCOOL))
    with pytest.raises(ValueError, match="needs a resolution"):
        mcool.fetch("chr4")
    with pytest.raises(ValueError, match="resolution 12345"):
        mcool.fetch("chr4", resolution=12345)


def test_errors_h5_and_hic():
    need(data(GSM_H5))
    h5 = hicx_matrix.open(data(GSM_H5))
    with pytest.raises(ValueError, match="available: none"):
        h5.fetch("chr1:0-1,000,000", normalization="KR")
    with pytest.raises(ValueError, match="obs_exp"):
        h5.fetch("chr1:0-1,000,000", transform="obs_exp")
    with pytest.raises(KeyError, match="unknown chromosome"):
        h5.fetch("chr2")
    hic = hicx_matrix.open(data(SMALL_HIC))
    with pytest.raises(ValueError, match="needs a resolution"):
        hic.fetch("NC_001133.9")
    with pytest.raises(ValueError, match="resolution 7"):
        hic.fetch("NC_001133.9", resolution=7)
    with pytest.raises(ValueError, match="unknown normalization 'FOO'"):
        hic.fetch("NC_001133.9", resolution=5000, normalization="FOO")
    with pytest.raises(KeyError, match="unknown chromosome 'All'"):
        hic.fetch("All", resolution=5000)
    with pytest.raises(FileNotFoundError):
        hicx_matrix.open(str(TEST_DATA / "does_not_exist.cool"))
    with pytest.raises(ValueError, match="neither"):
        hicx_matrix.open(__file__)
