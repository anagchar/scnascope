"""HiScanner (Zhao et al. 2025): single-cell copy number from scWGS read depth (and BAF).

Reads the output directory of the HiScanner pipeline (or its ``final_calls`` folder):

- ``final_calls/<cell>.txt``: one row per bin, ``CHROM START END RDR OBS EXP RDR_MEAN CN
  gamma``: read-depth ratio, observed and expected reads, the mean RDR of the bin's segment,
  the integer copy number and the cell's scaling factor (copy number ≈ RDR × gamma).
- ``final_calls/<cell>_seg.txt`` etc. are segment tables and are not needed.

Layers: ``total_cn`` (CN), ``rdr``, ``cn_scale`` (RDR × gamma: the bins on the copy-number
scale), ``read_counts`` (OBS); runs of the standard (RDR + BAF) pipeline also give ``baf``
(B / (A + B)) and ``baf_n`` (A + B reads), and its allele-specific calls (``CN`` "a|b",
``CN_A``/``CN_B``) ``major_cn``/``minor_cn``. Cells: ``gamma``. When cells were binned
differently, every cell is projected onto the most common bin grid.
"""

from __future__ import annotations

import re
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from ..genome import segments_to_bins
from ..model import CNVDataset, Modality, Semantics, chr_sort_key
from . import Reader, register

SKIP = re.compile(r"_(seg|seg_merged|seg_initial_anno|input_table)\.txt$|^all_")
HEAD = ("CHROM", "START", "END", "RDR", "CN")


def _calls_dir(path: Path) -> Path | None:
    for p in (path / "final_calls", path):
        if p.is_dir() and p.name == "final_calls":
            return p
    return None


def _cell_files(d: Path) -> list[Path]:
    return sorted(p for p in d.glob("*.txt") if not SKIP.search(p.name))


def _detect(path: Path) -> bool:
    d = _calls_dir(path)
    if d is None:
        return False
    for p in _cell_files(d)[:1]:
        with open(p) as f:
            head = f.readline().split()
        return all(h in head for h in HEAD)
    return False


def read_hiscanner(path: Path, genome: str = "hg38") -> CNVDataset:
    d = _calls_dir(path)
    if d is None:
        raise ValueError(f"no HiScanner final_calls folder under {path}")
    files = _cell_files(d)
    if not files:
        raise ValueError(f"no per-cell call files in {d}")

    # one small file per cell: read them in parallel (much faster on network drives)
    def read(p: Path) -> pd.DataFrame | None:
        t = pd.read_csv(p, sep="\t", dtype={"CHROM": str})
        if "CN_total" in t:  # standard pipeline: CN is "a|b", the total and alleles have columns
            t["CN"] = pd.to_numeric(t["CN_total"], errors="coerce")
        return t if all(h in t.columns for h in HEAD) else None

    with ThreadPoolExecutor(max_workers=8) as pool:  # more threads only slow network drives
        loaded = list(pool.map(read, files))
    tables = {
        p.name[: -len(".txt")]: t for p, t in zip(files, loaded, strict=True) if t is not None
    }
    cells = list(tables)

    grids = Counter(
        tuple(map(tuple, t[["CHROM", "START", "END"]].to_numpy().tolist())) for t in tables.values()
    )
    ref, n_ref = grids.most_common(1)[0]
    bins = pd.DataFrame(list(ref), columns=["chr", "start", "end"])
    bins = bins.assign(_k=bins["chr"].map(chr_sort_key)).sort_values(["_k", "start"])
    bins = bins.drop(columns="_k").reset_index(drop=True)
    key = pd.MultiIndex.from_frame(bins)

    cols = {"total_cn": "CN", "rdr": "RDR", "read_counts": "OBS"}
    # the standard (RDR + BAF) pipeline adds the phased allele reads A and B per bin
    allelic = all({"A", "B"} <= set(t.columns) for t in tables.values())
    if allelic:
        cols |= {"_a": "A", "_b": "B"}
    # ... and allele-specific calls CN_A / CN_B
    as_calls = all({"CN_A", "CN_B"} <= set(t.columns) for t in tables.values())
    if as_calls:
        cols |= {"_cna": "CN_A", "_cnb": "CN_B"}
    nb = len(bins)
    layers = {k: np.full((len(cells), nb), np.nan, dtype=np.float32) for k in cols}
    gamma = np.full(len(cells), np.nan)
    off_grid = []
    for i, c in enumerate(cells):
        t = tables[c]
        if "gamma" in t:
            gamma[i] = float(t["gamma"].iloc[0])
        idx = key.get_indexer(pd.MultiIndex.from_frame(t[["CHROM", "START", "END"]]))
        if (idx >= 0).all() and len(idx) == len(t):
            for k, src in cols.items():
                if src in t:
                    layers[k][i, idx] = t[src].to_numpy(np.float32)
        else:
            off_grid.append(c)
            seg = t.rename(columns={"CHROM": "chr", "START": "start", "END": "end"}).assign(cell=c)
            present = [src for src in cols.values() if src in t]
            proj = segments_to_bins(seg, bins, present, cells=[c])
            for k, src in cols.items():
                if src in proj:
                    layers[k][i] = proj[src][0]
    layers["cn_scale"] = (layers["rdr"] * gamma[:, None]).astype(np.float32)
    if as_calls:
        ca, cb = layers.pop("_cna"), layers.pop("_cnb")
        layers["major_cn"], layers["minor_cn"] = np.fmax(ca, cb), np.fmin(ca, cb)
    if allelic:  # BAF as B / (A + B), and the reads behind it
        a_, b_ = layers.pop("_a"), layers.pop("_b")
        with np.errstate(invalid="ignore", divide="ignore"):
            layers["baf"] = np.where(a_ + b_ > 0, b_ / (a_ + b_), np.nan).astype(np.float32)
        layers["baf_n"] = (a_ + b_).astype(np.float32)
    sem = {
        "total_cn": Semantics.INTEGER_CN,
        "rdr": Semantics.RATIO,
        "read_counts": Semantics.COUNT,
        "cn_scale": Semantics.CONTINUOUS_CN,
    }
    if allelic:
        sem |= {"baf": Semantics.BAF, "baf_n": Semantics.ALLELE_COUNT}
    if as_calls:
        sem |= {"major_cn": Semantics.ALLELE_CN, "minor_cn": Semantics.ALLELE_CN}
    obs = pd.DataFrame({"gamma": gamma}, index=cells)
    root = d.parent if d.name == "final_calls" else d
    return CNVDataset.build(
        cells=obs,
        bins=bins,
        layers=layers,
        semantics=sem,
        modality=Modality.SCDNA,
        tool="HiScanner",
        genome=genome,
        default_layer="total_cn",
        extra_uns={
            "source": root.name,
            "cells_projected_to_grid": off_grid,
            "cells_on_grid": n_ref,
        },
    )


def _read(path: Path, **opts) -> CNVDataset:
    return read_hiscanner(path, **opts)


register(Reader("hiscanner", "scDNA", _detect, _read, "HiScanner output directory (final_calls/)"))
