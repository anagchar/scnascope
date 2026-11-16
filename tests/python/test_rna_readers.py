"""CopyKAT and InferCNV outputs, in the layouts the real runs write."""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from scnascope.cli import main
from scnascope.model import Semantics
from scnascope.readers import detect, load

CELLS = ["AAAC", "GGTT", "TTAG", "CCGA"]


def write_copykat(d: Path, name: str = "test", python_port: bool = False) -> Path:
    """<name>_copykat_CNA_results.txt (windows x cells) and <name>_copykat_prediction.txt."""
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    # chr1: three windows, chr2: two; chromosome 23 is X
    for chrom, pos in ((1, 1042457), (1, 1262457), (1, 1482457), (2, 10000), (23, 300000)):
        vals = [0.1, -0.1, 0.0, 0.05] if chrom == 1 else [0.0, 0.0, 0.2, -0.2]
        rows.append(f"{chrom}\t{pos}\t{pos}\t" + "\t".join(map(str, vals)) + "\n")
    head = "chrom\tchrompos\tabspos\t" + "\t".join(CELLS) + "\n"
    (d / f"{name}_copykat_CNA_results.txt").write_text(head + rows[3] + "".join(rows[:3]) + rows[4])
    preds = ["aneuploid", "diploid", "aneuploid", "not.defined"]
    if python_port:
        preds = [
            "c2:aneuploid:high.conf",
            "c1:diploid:low.conf",
            "c2:aneuploid:high.conf",
            "c1:diploid:low.conf",
        ]
    (d / f"{name}_copykat_prediction.txt").write_text(
        "cell.names\tcopykat.pred\n"
        + "".join(f"{c}\t{p}\n" for c, p in zip(CELLS, preds, strict=True))
    )
    return d


def test_copykat(tmp_path):
    d = write_copykat(tmp_path / "CopyKat")
    assert detect(d).name == "copykat"
    ds = load(d)
    a = ds.adata
    assert a.uns["tool"] == "CopyKAT" and ds.semantics == {"expression": Semantics.CONTINUOUS_RATIO}
    assert a.var["chr"].tolist() == ["chr1"] * 3 + ["chr2", "chrX"]  # sorted, 23 -> X
    assert a.var["end"].tolist()[:2] == [1262457, 1482457]  # a window runs to the next one
    np.testing.assert_allclose(ds.layer("expression")[:, 0], [0.1, -0.1, 0.0, 0.05])
    assert a.obs["prediction"].tolist() == ["aneuploid", "diploid", "aneuploid", "not.defined"]
    assert "total" in ds.dendrograms


def test_copykat_python_port_prediction(tmp_path):
    ds = load(write_copykat(tmp_path / "res", "Thyroid", python_port=True))
    o = ds.adata.obs
    assert o["cluster"].tolist()[:2] == ["c2", "c1"]
    assert o["prediction"].tolist()[:2] == ["aneuploid", "diploid"]
    assert o["confidence"].tolist()[0] == "high.conf"


def test_copykat_several_runs_need_input_file(tmp_path):
    d = write_copykat(tmp_path / "ck", "a")
    write_copykat(d, "b")
    with pytest.raises(ValueError, match="several CopyKAT runs"):
        load(d, "copykat")
    assert load(d / "b_copykat_CNA_results.txt", "copykat").adata.uns["source"] == "b"


def write_infercnv(root: Path) -> Path:
    d = root / "output"
    d.mkdir(parents=True)
    (d / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings").write_text(
        "cell_group_name\tcell\n"
        "Thyrocytes.Thyrocytes_s1\tAAAC\nThyrocytes.Thyrocytes_s1\tGGTT\n"
        "Thyrocytes.Thyrocytes_s2\tTTAG\nT cell.T cell_s1\tCCGA\n"
    )
    head = "cell_group_name\tcnv_name\tstate\tchr\tstart\tend\n"
    # an earlier HMM step (not used) and the final, Bayesian-filtered one
    (d / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat").write_text(
        head + "Thyrocytes.Thyrocytes_s1\tchr1-region_1\t6\tchr1\t0\t248956422\n"
    )
    (
        d / "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"
    ).write_text(
        head
        + "Thyrocytes.Thyrocytes_s1\tchr1-region_1\t2\tchr1\t0\t50000000\n"
        + "Thyrocytes.Thyrocytes_s1\tchr2-region_1\t3\tchr2\t0\t10000000\n"
        + "Thyrocytes.Thyrocytes_s2\tchr2-region_2\t5\tchr2\t100000000\t120000000\n"
    )
    (d / "cell_annotations.txt").write_text(
        "AAAC\tThyrocytes\nGGTT\tThyrocytes\nTTAG\tThyrocytes\nCCGA\tT cell\nNOGR\tT cell\n"
    )
    return root


def test_infercnv(tmp_path):
    root = write_infercnv(tmp_path / "InferCNV")
    assert detect(root).name == "infercnv"
    ds = load(root)
    a = ds.adata
    assert a.uns["infercnv_regions"].startswith("HMM_CNV_predictions.")
    assert a.obs_names.tolist() == ["AAAC", "GGTT", "TTAG", "CCGA", "NOGR"]
    v = a.var
    T = ds.layer("total_cn")

    def w(c, mb):
        return int(
            np.flatnonzero((v["chr"] == c) & (v["start"] <= mb * 1e6) & (v["end"] > mb * 1e6))[0]
        )

    assert T[0, w("chr1", 10)] == 1 and T[1, w("chr1", 10)] == 1  # state 2 = one copy
    assert T[0, w("chr1", 100)] == 2  # outside the region: neutral
    assert T[2, w("chr2", 110)] == 4 and T[2, w("chr1", 10)] == 2  # s2: state 5 = 4 copies
    assert (T[3] == 2).all()  # a subcluster without regions
    assert np.isnan(T[4]).all()  # annotated but in no subcluster
    assert a.obs["cell_type"].tolist()[3] == "T cell"
    assert a.obs["subcluster"].tolist()[2] == "Thyrocytes.Thyrocytes_s2"
    assert a.obs["cnv_burden"].iloc[3] == 0 and a.obs["cnv_burden"].iloc[0] > 0


def test_cli_rna_tools(tmp_path):
    write_copykat(tmp_path / "ck")
    write_infercnv(tmp_path / "icnv")
    for tool, src in (("copykat", "ck"), ("infercnv", "icnv")):
        r = CliRunner().invoke(
            main,
            [
                "--mode",
                "rna",
                "--tool",
                tool,
                "-i",
                str(tmp_path / src),
                "-o",
                str(tmp_path / "o"),
                "--no-view",
            ],
        )
        assert r.exit_code == 0, r.output


@pytest.mark.skipif(not shutil.which("Rscript"), reason="needs R")
def test_infercnv_expression_from_the_r_object(tmp_path):
    """run.final.infercnv_obj, read as plain attributes (no infercnv package)."""
    root = write_infercnv(tmp_path / "InferCNV")
    obj = root / "output" / "run.final.infercnv_obj"
    subprocess.run(
        [
            "Rscript",
            "-e",
            "e <- matrix(c(1, 1.1, 0.5, 1, 1, 1, 2, 2), nrow = 2, "
            "dimnames = list(c('G1', 'G2'), c('AAAC', 'GGTT', 'TTAG', 'CCGA'))); "
            "go <- data.frame(chr = c('chr1', 'chr1'), start = c(100, 200), stop = c(150, 250), "
            "row.names = c('G1', 'G2')); "
            "x <- structure(list(), class = 'infercnv', expr.data = e, gene_order = go, "
            "reference_grouped_cell_indices = list(`T cell` = 4L), "
            "observation_grouped_cell_indices = list(Thyrocytes = 1:3)); "
            f"saveRDS(x, '{obj}')",
        ],
        check=True,
        capture_output=True,
    )
    ds = load(root)
    E = ds.layer("expression")
    assert ds.semantics["expression"] is Semantics.CONTINUOUS_RATIO
    # chr1's first window: the mean of log2 of the two genes
    np.testing.assert_allclose(E[:4, 0], [np.log2(1.1) / 2, -0.5, 0, 1], atol=1e-6)
    assert ds.adata.obs["reference"].tolist()[:4] == [False, False, False, True]
