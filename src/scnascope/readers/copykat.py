"""CopyKAT (Gao et al. 2021): aneuploid tumour cells from scRNA-seq copy-number ratios.

Reads a CopyKAT output folder (the R package or its Python port), or one of its
``<name>_copykat_CNA_results.txt`` files:

- ``<name>_copykat_CNA_results.txt``: ``chrom chrompos abspos <cells>``, one row per genomic
  window (220 kb by default; chromosome 23 is X), the cells' smoothed log copy-number ratios
- ``<name>_copykat_prediction.txt``: ``cell.names copykat.pred``: aneuploid / diploid (the
  Python port writes ``c1:diploid:low.conf``: cluster, prediction and confidence)

Layer ``expression`` (the ratios; windows past the heatmap's width are averaged in pairs). Cells:
``copykat`` (the prediction as written) and, when it has them, ``prediction``, ``cluster`` and
``confidence``. A folder with several runs (several ``*_CNA_results.txt``) is refused with
their names: choose one with ``--input-file``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..analysis import ward_manhattan
from ..genome import coarsen_bins
from ..model import CNVDataset, Modality, Semantics, chr_sort_key
from . import Reader, register

SUFFIX = "_copykat_CNA_results.txt"
#: Ward on thousands of windows is slow: the cells are clustered on blocks of this many columns
CLUSTER_COLS = 600


def _results(path: Path) -> list[Path]:
    if path.is_file():
        return [path] if path.name.endswith(SUFFIX) else []
    if path.is_dir():
        return sorted(p for p in path.glob(f"*{SUFFIX}") if not p.name.startswith("._"))
    return []


def _detect(path: Path) -> bool:
    found = _results(path)
    if not found:
        return False
    with open(found[0]) as f:
        return f.readline().split("\t")[:3] == ["chrom", "chrompos", "abspos"]


def _chr(c: str) -> str:
    c = str(c).removeprefix("chr")
    return "chr" + {"23": "X", "24": "Y"}.get(c, c)


def read_copykat(path: Path, genome: str = "hg38") -> CNVDataset:
    found = _results(path)
    if not found:
        raise ValueError(f"no CopyKAT *{SUFFIX} in {path}")
    if len(found) > 1:
        names = ", ".join(p.name for p in found)
        raise ValueError(
            f"{path} holds several CopyKAT runs ({names}); choose one with --input-file"
        )
    res = found[0]
    with open(res) as f:
        head = f.readline().rstrip("\n").split("\t")
    cells = head[3:]
    types = {c: np.float32 for c in cells} | {"chrom": str}
    df = pd.read_csv(res, sep="\t", dtype=types, engine="c")
    chrs = df["chrom"].map(_chr).to_numpy()
    pos = df["chrompos"].to_numpy(np.int64)
    rank = {c: i for i, c in enumerate(sorted(set(chrs), key=chr_sort_key))}
    order = np.lexsort((pos, np.array([rank[c] for c in chrs])))
    chrs, pos = chrs[order], pos[order]
    X = df[cells].to_numpy(np.float32)[order].T  # cells x windows
    del df
    # windows are given by one position: each runs to the next window of its chromosome
    step = int(np.median(np.diff(pos)[np.diff(pos) > 0])) if len(pos) > 1 else 220_000
    end = np.empty_like(pos)
    for i in range(len(pos)):
        last = i + 1 == len(pos) or chrs[i + 1] != chrs[i]
        end[i] = pos[i] + step if last else pos[i + 1]
    bins = pd.DataFrame({"chr": chrs, "start": pos, "end": end})
    bins, X = coarsen_bins(bins, X)

    obs = pd.DataFrame(index=cells)
    pred = res.with_name(res.name[: -len(SUFFIX)] + "_copykat_prediction.txt")
    if pred.is_file():
        p = pd.read_csv(pred, sep="\t", dtype=str).drop_duplicates("cell.names")
        p = p.set_index("cell.names")["copykat.pred"].reindex(cells)
        obs["copykat"] = pd.Categorical(p.to_numpy(dtype=object))
        parts = p.str.split(":", expand=True)
        if parts.shape[1] == 3:  # the Python port: cluster:prediction:confidence
            for k, name in enumerate(("cluster", "prediction", "confidence")):
                obs[name] = pd.Categorical(parts[k].to_numpy(dtype=object))
        else:
            obs["prediction"] = obs["copykat"]

    dendros = {}
    if 3 <= len(cells) <= 20000:
        d = ward_manhattan(X, max_cols=CLUSTER_COLS)
        method = f"ward.D2, Manhattan on the log ratios (averaged to {CLUSTER_COLS} columns)"
        dendros["total"] = d | {"layers": ["expression"], "method": method}
    run = res.name[: -len(SUFFIX)]
    return CNVDataset.build(
        cells=obs,
        bins=bins,
        layers={"expression": X},
        semantics={"expression": Semantics.CONTINUOUS_RATIO},
        modality=Modality.SCRNA,
        tool="CopyKAT",
        genome=genome,
        default_layer="expression",
        extra_uns={
            "source": run or res.parent.name,
            "dendrograms": dendros,
            # expression-based ratios: the DNA QC metrics do not apply
            "qc": {"version": 1, "computed": []},
        },
    )


def _read(path: Path, **opts) -> CNVDataset:
    return read_copykat(path, **opts)


register(Reader("copykat", "scRNA", _detect, _read, "CopyKAT output (*_copykat_CNA_results.txt)"))
