"""CHISEL, HiScanner and ASCAT.sc-from-R readers, on small files in each tool's own format
(column layouts copied from real outputs)."""

import shutil
from pathlib import Path

import numpy as np
import pytest

from scnascope.model import Semantics
from scnascope.prepare import prepare
from scnascope.readers import detect, load

CHISEL_HEAD = (
    "#CHR\tSTART\tEND\tCELL\tNORM_COUNT\tCOUNT\tRDR\tA_COUNT\tB_COUNT\tBAF\tCLUSTER\t"
    "HAP_CN\tCORRECTED_HAP_CN\n"
)


def write_chisel(root: Path) -> None:
    (root / "calls").mkdir(parents=True)
    (root / "clones").mkdir()
    rows = []
    for cell, haps in (
        ("AAAC", ["1|1", "2|1", "1|0"]),
        ("GGTT", ["1|1", "1|1", "1|1"]),
        ("TTAG", ["2|2", "1|2", "0|0"]),
    ):
        for i, (chrom, s) in enumerate((("chr1", 0), ("chr1", 5000000), ("chr2", 0))):
            a, b = map(int, haps[i].split("|"))
            rdr = (a + b) / 2
            rows.append(
                f"{chrom}\t{s}\t{s + 5000000}\t{cell}\t1000\t{int(100 * rdr)}\t{rdr}\t5\t5\t0.5"
                f"\t17\t{haps[i]}\t{haps[i]}\n"
            )
    # rows of chr2 first, to check bins are sorted
    (root / "calls" / "calls.tsv").write_text(
        CHISEL_HEAD + "".join(rows[2::3] + [r for k, r in enumerate(rows) if k % 3 != 2])
    )
    (root / "clones" / "mapping.tsv").write_text(
        "#CELL\tCLUSTER\tCLONE\nAAAC\t17\tClone1\nGGTT\t17\tNone\nTTAG\t17\tClone1\n"
    )


def test_chisel(tmp_path):
    write_chisel(tmp_path)
    assert detect(tmp_path).name == "chisel"
    ds = prepare(load(tmp_path))
    a = ds.adata
    assert a.obs.index.tolist() == ["AAAC", "GGTT", "TTAG"]
    assert a.var[["chr", "start"]].values.tolist() == [["chr1", 0], ["chr1", 5000000], ["chr2", 0]]
    np.testing.assert_array_equal(ds.layer("total_cn")[0], [2, 3, 1])
    np.testing.assert_array_equal(ds.layer("major_cn")[2], [2, 2, 0])
    np.testing.assert_array_equal(ds.layer("minor_cn")[2], [2, 1, 0])
    np.testing.assert_array_equal(ds.layer("hap_a_cn")[2], [2, 1, 0])
    assert ds.semantics["baf"] is Semantics.BAF and ds.semantics["rdr"] is Semantics.RATIO
    # the reads behind the BAF: A_COUNT + B_COUNT (5 + 5 in every bin)
    assert ds.semantics["baf_n"] is Semantics.ALLELE_COUNT and (ds.layer("baf_n") == 10).all()
    assert a.obs["clone"].tolist()[0] == "Clone1" and a.obs["clone"].isna().tolist()[1]
    assert "allele" in ds.dendrograms


def write_hiscanner(root: Path) -> None:
    d = root / "final_calls"
    d.mkdir(parents=True)
    head = "CHROM\tSTART\tEND\tRDR\tOBS\tEXP\tRDR_MEAN\tCN\tgamma\n"
    grid = [("chr1", 10033, 1975218), ("chr1", 1975219, 3115514), ("chr2", 10005, 2000000)]
    for cell, cns, gamma in (("SRR1", [2, 3, 2], 2.0), ("SRR2", [2, 2, 1], 2.0)):
        rows = [
            f"{c}\t{s}\t{e}\t{cn / gamma}\t100\t100\t1\t{cn}.0\t{gamma}\n"
            for (c, s, e), cn in zip(grid, cns, strict=True)
        ]
        (d / f"{cell}.txt").write_text(head + "".join(rows))
        (d / f"{cell}_seg.txt").write_text("chrom\tstart\tend\n")
    # a cell binned differently: one bin split in two
    rows = [
        "chr1\t10033\t1975218\t1.0\t100\t100\t1\t2.0\t2.0\n",
        "chr1\t1975219\t2500000\t2.0\t100\t100\t1\t4.0\t2.0\n",
        "chr1\t2500001\t3115514\t2.0\t100\t100\t1\t4.0\t2.0\n",
        "chr2\t10005\t2000000\t1.0\t100\t100\t1\t2.0\t2.0\n",
    ]
    (d / "SRR3.txt").write_text(head + "".join(rows))


def test_hiscanner(tmp_path):
    write_hiscanner(tmp_path)
    assert detect(tmp_path).name == "hiscanner"
    ds = prepare(load(tmp_path))
    a = ds.adata
    assert a.obs.index.tolist() == ["SRR1", "SRR2", "SRR3"]
    assert a.n_vars == 3
    np.testing.assert_array_equal(ds.layer("total_cn")[0], [2, 3, 2])
    np.testing.assert_array_equal(ds.layer("total_cn")[2], [2, 4, 2])  # projected by midpoint
    np.testing.assert_allclose(ds.layer("cn_scale")[0], [2, 3, 2])  # RDR × gamma
    assert a.uns["cells_projected_to_grid"] == ["SRR3"]


def test_hiscanner_standard_pipeline_baf(tmp_path):
    """The standard (RDR + BAF) pipeline: phased allele reads A, B and allele-specific calls
    (CN written "a|b", with CN_A / CN_B / CN_total)."""
    d = tmp_path / "final_calls"
    d.mkdir()
    head = (
        "CHROM\tSTART\tEND\tA\tB\tN\tTOTAL\tpBAF\tBAF\tRDR\tOBS\tEXP\tRDR_MEAN"
        "\tCN\tCN_A\tCN_B\tCN_total\tgamma\n"
    )
    rows = [
        "chr1\t1\t1000000\t3\t1\t2\t4\t0.25\t0.25\t1\t50\t50\t1\t2|1\t2.0\t1.0\t3.0\t2.0\n",
        "chr1\t1000001\t2000000\t0\t0\t0\t0\t0\t0\t1\t50\t50\t1\t1|1\t1.0\t1.0\t2.0\t2.0\n",
    ]
    for cell in ("C1", "C2"):
        (d / f"{cell}.txt").write_text(head + "".join(rows))
    (d / "all_cells_bins.txt").write_text(head + "".join(rows))  # a summary, not a cell
    ds = load(tmp_path)
    assert ds.adata.obs_names.tolist() == ["C1", "C2"]
    np.testing.assert_allclose(ds.layer("baf")[0, 0], 0.25)  # B / (A + B)
    assert np.isnan(ds.layer("baf")[0, 1]) and ds.layer("baf_n")[0].tolist() == [4, 0]
    assert ds.layer("total_cn")[0].tolist() == [3, 2]  # CN_total, not the "a|b" string
    assert ds.layer("major_cn")[0].tolist() == [2, 1] and ds.layer("minor_cn")[0].tolist() == [1, 1]
    assert ds.semantics["baf_n"] is Semantics.ALLELE_COUNT and "_a" not in ds.layer_names


def test_hiscanner_rdr_only_has_no_baf(tmp_path):
    write_hiscanner(tmp_path)
    assert "baf" not in load(tmp_path).layer_names


@pytest.mark.skipif(not shutil.which("Rscript"), reason="needs R")
def test_ascatsc_rds_failure_is_reported(tmp_path):
    p = tmp_path / "x.rds"
    p.write_bytes(b"not an rds")
    assert detect(p).name == "ascatsc_r"
    with pytest.raises(RuntimeError, match="R export of x.rds failed"):
        load(p)
