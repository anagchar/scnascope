"""InferCNV (Trinity CTAT): copy-number states of tumour subclusters from scRNA-seq expression.

Reads an InferCNV ``out_dir`` (or a folder holding ``output/``) through its plain-text HMM
results, so no R is needed:

- ``*.pred_cnv_regions.dat``: ``cell_group_name cnv_name state chr start end``, the HMM state
  of each region of each subcluster; the Bayesian-filtered ``HMM_CNV_predictions.*`` is used
  when present, otherwise the last ``NN_HMM_pred*`` file
- ``*.cell_groupings``: ``cell_group_name cell``, the cells of each subcluster
- ``cell_annotations.txt`` (optional): the cell types the run was given

HMM i6 states as copy number: 1 → 0 (complete loss), 2 → 1, 3 → 2 (neutral), 4 → 3, 5 → 4,
6 → 6 (> 2x). Each cell is painted with its subcluster's regions on 1 Mb windows over a
neutral background (``total_cn``). Cells: ``cell_type``, ``subcluster``, ``reference``.

With R on the PATH, the smoothed expression InferCNV plots (``infercnv.png``) is read too, from
``run.final.infercnv_obj`` (R/export_infercnv.R: the object's attributes, so the infercnv
package is not needed): ``expression``, its genes averaged per window, as log2 of InferCNV's
values (1 = the references' level), and the reference cells.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from ..analysis import ward_manhattan
from ..genome import CHROM_SIZES
from ..model import CNVDataset, Modality, Semantics, chr_sort_key
from . import Reader, register

WINDOW = 1_000_000
SCRIPT = Path(__file__).parents[1] / "R" / "export_infercnv.R"
#: the objects with InferCNV's final smoothed expression, best first
OBJECTS = ("run.final.infercnv_obj", "preliminary.infercnv_obj")
STATE_CN = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 6}


def _out_dir(path: Path) -> Path | None:
    for d in (path, path / "output"):
        if d.is_dir() and any(d.glob("*.pred_cnv_regions.dat")) and any(d.glob("*.cell_groupings")):
            return d
    return None


def _detect(path: Path) -> bool:
    return _out_dir(path) is not None


def _step(p: Path) -> int:
    m = re.match(r"(\d+)_", p.name)
    return int(m[1]) if m else -1


def _files(d: Path) -> tuple[Path, Path]:
    regions = sorted(d.glob("*.pred_cnv_regions.dat"))
    final = [p for p in regions if p.name.startswith("HMM_CNV_predictions.")]
    reg = final[0] if final else max(regions, key=_step)
    groupings = sorted(d.glob("*.cell_groupings"), key=_step)
    return reg, groupings[-1]


def _expression(
    d: Path, cells: list[str], bins: pd.DataFrame, rscript: str | None = None
) -> tuple[np.ndarray, pd.DataFrame] | None:
    """InferCNV's smoothed expression (cells x windows, log2) and its groups, or None when
    there is no object or no R."""
    obj = next((d / n for n in OBJECTS if (d / n).is_file()), None)
    exe = rscript or shutil.which("Rscript")
    if obj is None or not exe:
        return None
    print(f"scnascope: reading {obj.name} with R ...", file=sys.stderr)
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run([exe, str(SCRIPT), str(obj), tmp], capture_output=True, text=True)
        out = Path(tmp)
        if r.returncode != 0 or not (out / "expr.f32").is_file():
            tail = "\n".join((r.stderr or r.stdout).strip().splitlines()[-5:])
            print(f"scnascope: InferCNV expression not read: {tail}", file=sys.stderr)
            return None
        genes = pd.read_csv(out / "genes.tsv", sep="\t", dtype={"gene": str, "chr": str})
        ecells = (out / "cells.txt").read_text().splitlines()
        E = np.fromfile(out / "expr.f32", dtype="<f4").reshape(len(ecells), len(genes))
        groups = pd.read_csv(out / "groups.tsv", sep="\t", dtype=str)
    # genes onto the windows by their start
    key = pd.MultiIndex.from_arrays([bins["chr"].to_numpy(), (bins["start"] // WINDOW).to_numpy()])
    gk = pd.MultiIndex.from_arrays(
        [genes["chr"].to_numpy(), (genes["start"].astype(np.int64) // WINDOW).to_numpy()]
    )
    col = key.get_indexer(gk)
    ok = col >= 0
    with np.errstate(divide="ignore", invalid="ignore"):
        V = np.log2(E[:, ok]).astype(np.float64)
    S = np.zeros((len(ecells), len(bins)))
    C = np.zeros((len(ecells), len(bins)))
    fin = np.isfinite(V)
    np.add.at(S.T, col[ok], np.where(fin, V, 0).T)
    np.add.at(C.T, col[ok], fin.T)
    with np.errstate(invalid="ignore", divide="ignore"):
        W = np.where(C > 0, S / C, np.nan)
    rows = pd.Index(ecells).get_indexer(cells)
    X = np.full((len(cells), len(bins)), np.nan, dtype=np.float32)
    X[rows >= 0] = W[rows[rows >= 0]]
    return X, groups.drop_duplicates("cell").set_index("cell")


def read_infercnv(path: Path, genome: str = "hg38", rscript: str | None = None) -> CNVDataset:
    d = _out_dir(path)
    if d is None:
        raise ValueError(
            f"no InferCNV HMM results (*.pred_cnv_regions.dat, *.cell_groupings) in {path}"
        )
    reg_file, grp_file = _files(d)
    reg = pd.read_csv(reg_file, sep="\t", dtype={"cell_group_name": str, "chr": str})
    grp = pd.read_csv(grp_file, sep="\t", dtype=str).drop_duplicates("cell")
    ann_file = next(
        (f for f in (d / "cell_annotations.txt", d.parent / "cell_annotations.txt") if f.is_file()),
        None,
    )
    ann = (
        pd.read_csv(ann_file, sep="\t", header=None, names=["cell", "cell_type"], dtype=str)
        .drop_duplicates("cell")
        .set_index("cell")["cell_type"]
        if ann_file
        else None
    )
    cells = grp["cell"].tolist()
    if ann is not None:
        cells += [c for c in ann.index if c not in set(cells)]

    # 1 Mb windows over the chromosomes the run has (the genome's lengths, or the regions' ends)
    sizes = CHROM_SIZES.get(genome, {})
    chroms = sorted({str(c) for c in reg["chr"]}, key=chr_sort_key)
    rows = []
    for c in chroms:
        n = (
            sizes.get(c)
            or sizes.get(c.removeprefix("chr"))
            or int(reg.loc[reg["chr"] == c, "end"].max())
        )
        rows += [(c, k * WINDOW, min((k + 1) * WINDOW, n)) for k in range(int(np.ceil(n / WINDOW)))]
    bins = pd.DataFrame(rows, columns=["chr", "start", "end"])
    mid = ((bins["start"] + bins["end"]) / 2).to_numpy()
    chr_b = bins["chr"].to_numpy()

    # each subcluster's states, painted on a neutral background
    states = {}
    for g, r in reg.groupby("cell_group_name"):
        s = np.full(len(bins), 3.0)
        for c, a, b, st in r[["chr", "start", "end", "state"]].itertuples(index=False):
            s[(chr_b == c) & (mid >= a) & (mid <= b)] = st
        states[g] = s
    group_of = grp.set_index("cell")["cell_group_name"]
    neutral = np.full(len(bins), 3.0)
    S = np.vstack(
        [
            states.get(group_of.get(c), neutral)
            if c in group_of.index
            else np.full(len(bins), np.nan)
            for c in cells
        ]
    )
    cn = np.vectorize(
        lambda v: STATE_CN.get(int(v), np.nan) if np.isfinite(v) else np.nan, otypes=[float]
    )(S)

    obs = pd.DataFrame(index=cells)
    if ann is not None:
        obs["cell_type"] = pd.Categorical(ann.reindex(cells).to_numpy(dtype=object))
    obs["subcluster"] = pd.Categorical(group_of.reindex(cells).to_numpy(dtype=object))
    with np.errstate(invalid="ignore"):
        obs["cnv_burden"] = np.where(np.isfinite(S).any(axis=1), np.nanmean(S != 3, axis=1), np.nan)

    layers = {"total_cn": cn.astype(np.float32)}
    sem = {"total_cn": Semantics.INTEGER_CN}
    expr = _expression(d, cells, bins, rscript)
    if expr is not None:
        layers["expression"], groups = expr
        sem["expression"] = Semantics.CONTINUOUS_RATIO
        kind = groups["kind"].reindex(cells)
        obs["reference"] = kind.eq("reference").to_numpy()
        if "cell_type" not in obs:
            obs["cell_type"] = pd.Categorical(groups["group"].reindex(cells).to_numpy(dtype=object))

    dendros = {}
    if 3 <= len(cells) <= 20000:
        dd = ward_manhattan(cn, max_cols=600)
        dendros["total"] = dd | {
            "layers": ["total_cn"],
            "method": "ward.D2, Manhattan on the HMM states",
        }
    return CNVDataset.build(
        cells=obs,
        bins=bins,
        layers=layers,
        semantics=sem,
        modality=Modality.SCRNA,
        tool="InferCNV",
        genome=genome,
        default_layer="total_cn",
        extra_uns={
            "source": (d.parent.name if d.name == "output" else d.name),
            "infercnv_regions": reg_file.name,
            "infercnv_expression": "run.final.infercnv_obj" if expr is not None else "",
            "dendrograms": dendros,
            "qc": {"version": 1, "computed": ["cnv_burden"]},
        },
    )


def _read(path: Path, **opts) -> CNVDataset:
    return read_infercnv(path, **opts)


register(
    Reader("infercnv", "scRNA", _detect, _read, "InferCNV out_dir (HMM *.pred_cnv_regions.dat)")
)
