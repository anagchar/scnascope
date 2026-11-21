"""epiAneufinder's results_table.tsv, cell annotations, and the --mode atac run."""

import shutil
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from scnascope.annotate import annotate
from scnascope.cli import main
from scnascope.model import Semantics
from scnascope.readers import detect, load

CELLS = ["AAACGAAAGC-1", "TTGCAGGTCA-1", "GGGTCATTAC-1"]
# chr1 windows with one filtered out (3-4 Mb), and chr2
WINDOWS = [("chr1", 1, 1e6), ("chr1", 1e6 + 1, 2e6), ("chr1", 2e6 + 1, 3e6), ("chr1", 4e6 + 1, 5e6)]
WINDOWS += [("chr2", 1, 1e6), ("chr2", 1e6 + 1, 2e6)]
CALLS = np.array(  # windows x cells, 0 loss / 1 normal / 2 gain
    [[1, 2, 0], [1, 2, 0], [1, 2, 1], [1, 1, 1], [1, 1, 0], [1, 1, 0]]
)


def r_num(v: float) -> str:
    # R's write.table prints 1e+06 for one million
    return f"{v:.0e}" if v >= 1e5 and v % 1e5 == 0 else str(int(v))


def write_epi(root: Path) -> Path:
    """An epiAneufinder output folder, written like write.table(..., quote = FALSE)."""
    d = root / "epiAneufinder_results"
    d.mkdir(parents=True)
    head = "seq start end " + " ".join(f"cell-{c}" for c in CELLS)
    rows = [
        f"{i + 1} {c} {r_num(s)} {r_num(e)} " + " ".join(map(str, CALLS[i]))
        for i, (c, s, e) in enumerate(WINDOWS)
    ]
    (d / "results_table.tsv").write_text("\n".join([head, *rows]) + "\n")
    return root


def test_reads_results_table(tmp_path):
    root = write_epi(tmp_path / "epi_out")
    assert detect(root).name == "epianeufinder"
    ds = load(root)
    a = ds.adata
    assert a.uns["tool"] == "epiAneufinder" and a.uns["modality"] == "scATAC"
    assert a.obs_names.tolist() == CELLS
    assert ds.semantics == {"state": Semantics.TERNARY_STATE}
    assert a.var["end"].tolist()[:2] == [1_000_000, 2_000_000]
    np.testing.assert_array_equal(ds.layer("state"), CALLS.T - 1)
    np.testing.assert_allclose(a.obs["cnv_burden"], [0, 3 / 6, 4 / 6])
    np.testing.assert_allclose(a.obs["gain_fraction"], [0, 3 / 6, 0])
    np.testing.assert_allclose(a.obs["loss_fraction"], [0, 0, 4 / 6])


def test_reads_the_table_file_and_results_folder(tmp_path):
    root = write_epi(tmp_path / "epi_out")
    for p in (root / "epiAneufinder_results", root / "epiAneufinder_results" / "results_table.tsv"):
        assert load(p, "epianeufinder").adata.n_obs == 3


@pytest.mark.parametrize(
    "barcodes",
    [CELLS, [f"cell-{c}" for c in CELLS], [f"sample1#{c[:-2]}" for c in CELLS]],
)
def test_annotations_match_barcodes(tmp_path, barcodes):
    ds = load(write_epi(tmp_path / "epi"))
    f = tmp_path / "annot.csv"
    pd.DataFrame({"cell": barcodes, "annot": ["normal", "tumour", "tumour"]}).to_csv(f, index=False)
    assert annotate(ds, f) == 3
    assert ds.adata.obs["cell_type"].tolist() == ["normal", "tumour", "tumour"]


def test_annotations_that_match_nothing_are_an_error(tmp_path):
    ds = load(write_epi(tmp_path / "epi"))
    f = tmp_path / "annot.tsv"
    f.write_text("barcode\tcelltype\nXXXX\tB cell\n")
    with pytest.raises(ValueError, match="no barcode"):
        annotate(ds, f)


def test_atac_run_with_annotations(tmp_path):
    root = write_epi(tmp_path / "epi_out")
    f = tmp_path / "cells.tsv"
    f.write_text(
        "barcode\tcell_type\tsample\n" + "".join(f"{c}\tT{i}\tS1\n" for i, c in enumerate(CELLS))
    )
    out = tmp_path / "results"
    r = CliRunner().invoke(
        main,
        [
            "--mode",
            "atac",
            "--epianeufinder",
            "--input",
            str(root),
            "--annotations",
            str(f),
            "--outdir",
            str(out),
            "--no-view",
        ],
    )
    assert r.exit_code == 0, r.output
    assert "annotations for 3 of 3 cells" in r.output
    a = ad.read_h5ad(out / "epi_out.h5ad")
    assert a.obs["cell_type"].tolist() == ["T0", "T1", "T2"] and a.uns["tool"] == "epiAneufinder"
    cells = pd.read_csv(out / "epi_out_cells.tsv", sep="\t")
    assert {"cnv_burden", "cell_type", "sample"} <= set(cells.columns)


def test_atac_tool_needs_atac_mode(tmp_path):
    root = write_epi(tmp_path / "epi_out")
    r = CliRunner().invoke(
        main,
        ["--mode", "dna", "--tool", "epianeufinder", "-i", str(root), "-o", str(tmp_path / "o")],
    )
    assert r.exit_code != 0 and "scATAC tool" in r.output


@pytest.mark.skipif(not shutil.which("Rscript"), reason="needs R")
def test_gc_corrected_counts_from_the_rds(tmp_path):
    """counts_gc_corrected.rds (windows x cells, a data.table) becomes reads per bin."""
    root = write_epi(tmp_path / "epi_out")
    rds = root / "epiAneufinder_results" / "counts_gc_corrected.rds"
    # columns in another order than the table's, to check they are matched by name
    cols = ", ".join(
        f"`cell-{c}` = c({', '.join(str(10 * k + w) for w in range(len(WINDOWS)))})"
        for k, c in reversed(list(enumerate(CELLS)))
    )
    subprocess.run(
        ["Rscript", "-e", f"x <- data.frame({cols}, check.names = FALSE); saveRDS(x, '{rds}')"],
        check=True,
        capture_output=True,
    )
    ds = load(root)
    assert ds.semantics["read_counts"] is Semantics.COUNT
    R = ds.layer("read_counts")
    np.testing.assert_array_equal(R[1], [10, 11, 12, 13, 14, 15])  # the second cell, by name


def test_bin_size_merges_windows(tmp_path):
    """1 Mb on 1 Mb windows is refused; 2 Mb merges chr1's 0-2 Mb, 2-3(+4-5) Mb and chr2."""
    root = write_epi(tmp_path / "epi_out")
    with pytest.raises(ValueError, match="not larger"):
        load(root, "epianeufinder", bin_size="1Mb")
    ds = load(root, "epianeufinder", bin_size="2Mb")
    v = ds.adata.var
    # chr1 windows start at 1, 1e6+1, 2e6+1, 4e6+1, chr2 at 1, 1e6+1: 2 Mb cells 0 0 1 2 | 0 0
    assert v["chr"].tolist() == ["chr1", "chr1", "chr1", "chr2"]
    assert v["n_windows"].tolist() == [2, 1, 1, 2]
    assert v["start"].tolist()[:2] == [1, 2_000_001] and v["end"].tolist()[0] == 2_000_000
    S = ds.layer("state")
    # cell 2: gain, gain in the first merged bin; cell 3: loss, loss then normal; chr2: loss, loss
    np.testing.assert_array_equal(S[:, 0], [0, 1, -1])
    np.testing.assert_array_equal(S[:, 3], [0, 0, -1])
    # the cells' burden is still epiAneufinder's own windows'
    np.testing.assert_allclose(ds.adata.obs["cnv_burden"], [0, 3 / 6, 4 / 6])
    assert ds.adata.uns["bin_size"] == 2_000_000


def test_bin_size_ties_are_normal():
    from scnascope.readers.epianeufinder import merge_windows, parse_bin_size

    assert parse_bin_size("1Mb") == 1_000_000 and parse_bin_size("500kb") == 500_000
    assert parse_bin_size("2e6") == 2_000_000
    bins = pd.DataFrame(
        {"chr": ["chr1"] * 4, "start": [1, 101, 201, 301], "end": [100, 200, 300, 400]}
    )
    state = np.array([[1, -1, 0, 0], [1, 1, -1, 0], [1, -1, np.nan, np.nan]], dtype=np.float32)
    counts = np.array([[1, 2, 3, 4], [0, 0, 0, 0], [5, np.nan, 1, 1]], dtype=np.float32)
    b, S, C = merge_windows(bins, state, counts, 200)
    assert b["n_windows"].tolist() == [2, 2]
    # ties (gain/loss, normal/loss) are normal; two gains are a gain; all missing is NaN
    np.testing.assert_array_equal(S[:, 0], [0, 1, 0])
    assert S[1, 1] == 0 and np.isnan(S[2, 1])
    np.testing.assert_array_equal(C, [[3, 7], [0, 0], [5, 2]])


def test_bin_size_cli_reruns_for_a_new_size(tmp_path):
    root = write_epi(tmp_path / "epi_out")
    out = tmp_path / "o"
    base = ["--mode", "atac", "--epianeufinder", "-i", str(root), "-o", str(out), "--no-view"]
    assert CliRunner().invoke(main, base).exit_code == 0
    r = CliRunner().invoke(main, [*base, "--bin-size", "2Mb"])
    assert r.exit_code == 0, r.output
    assert "unchanged" not in r.output and ad.read_h5ad(out / "epi_out.h5ad").n_vars == 4
    r = CliRunner().invoke(main, [*base, "--bin-size", "2Mb"])
    assert "unchanged" in r.output
    r = CliRunner().invoke(
        main, ["--mode", "dna", "--chisel", "-i", str(root), "-o", str(out), "--bin-size", "1Mb"]
    )
    assert r.exit_code != 0 and "epiAneufinder runs only" in r.output


def test_logr_against_the_cells_own_type():
    """log2(share of the cell's fragments in a bin / median share in its cell type), centred."""
    from scnascope.model import CNVDataset, Modality
    from scnascope.readers.epianeufinder import add_logr

    # two cell types with different accessibility in bin 2; cell 3 has lost bin 3
    types = ["A"] * 3 + ["B"] * 3
    C = np.array(
        [[10, 10, 40, 10, 10]] * 3
        + [[10, 10, 10, 10, 10], [10, 10, 10, 10, 10], [10, 10, 10, 0, 10]],
        dtype=np.float32,
    )
    bins = pd.DataFrame({"chr": ["chr1"] * 5, "start": np.arange(5) * 100 + 1})
    bins["end"] = bins["start"] + 99
    ds = CNVDataset.build(
        cells=pd.DataFrame({"cell_type": types}, index=[f"c{i}" for i in range(6)]),
        bins=bins,
        layers={"state": np.zeros((6, 5), np.float32), "read_counts": C},
        semantics={"state": Semantics.TERNARY_STATE, "read_counts": Semantics.COUNT},
        modality=Modality.SCATAC,
        tool="epiAneufinder",
        genome="hg38",
        default_layer="state",
    )
    ref = add_logr(ds, min_cells=3)
    assert ref.startswith("its own cell_type (2 with >= 3 cells")
    L = ds.layer("logr")
    # type A's open bin 2 is its own normal: ~0, not a gain
    assert abs(L[0, 2]) < 0.05
    assert L[5, 3] < -3  # the lost bin
    assert ds.semantics["logr"] is Semantics.CONTINUOUS_RATIO and "mapd" in ds.adata.obs
    # against all cells instead, type A's open bin looks like a gain
    add_logr(ds, by="none")
    assert ds.layer("logr")[0, 2] > 0.5
