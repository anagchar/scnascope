"""epiAneufinder (Ramakrishnan et al. 2023): gains and losses from scATAC-seq read depth.

Reads the run's output folder (``epiAneufinder_results/``, or the folder above it):

- ``results_table.tsv``: R ``write.table`` output (space-separated, with row names), columns
  ``seq start end cell-<barcode> ...``, one row per genomic window, values 0 = loss, 1 =
  normal, 2 = gain.

Layer ``state`` (-1 loss, 0 normal, +1 gain). Cells (barcodes without the ``cell-`` prefix):
``cnv_burden`` (fraction of windows not normal, epiAneufinder's ``cnv_burden_per_cell``),
``gain_fraction`` and ``loss_fraction``.

With R on the PATH, ``counts_gc_corrected.rds`` (next to the table: the GC-corrected fragments
per window that epiAneufinder segments, same windows and order) is read too, through
R/export_epianeufinder.R: layer ``read_counts``, shown as reads per bin.

``bin_size`` (e.g. 1 Mb) merges the windows into larger bins on a grid of that size, per
chromosome: fragment counts are summed, and each merged bin takes the call most of its windows
have (a tie with normal, or between gain and loss, is normal). A merged bin spans its first
window's start to its last window's end, and ``var["n_windows"]`` says how many windows it holds
(fewer at chromosome ends and where epiAneufinder dropped windows): no window is invented. The
cells' ``cnv_burden`` and gain/loss fractions stay those of epiAneufinder's own windows.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from ..analysis import ward_d_manhattan
from ..model import CNVDataset, Modality, Semantics, chr_sort_key
from . import Reader, register

TABLE = "results_table.tsv"
PREFIX = "cell-"
COUNTS = "counts_gc_corrected.rds"
SCRIPT = Path(__file__).parents[1] / "R" / "export_epianeufinder.R"


def _counts(rds: Path, cols: list[str], n_win: int, rscript: str | None) -> np.ndarray | None:
    """GC-corrected fragment counts (cells x windows, in the table's window order) for the
    table's cells `cols`, or None when there is no R or the file does not match."""
    exe = rscript or shutil.which("Rscript")
    if not exe:
        print(f"scnascope: no Rscript; {rds.name} (reads per bin) not read", file=sys.stderr)
        return None
    print(f"scnascope: reading {rds.name} with R ...", file=sys.stderr)
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run([exe, str(SCRIPT), str(rds), tmp], capture_output=True, text=True)
        out = Path(tmp)
        if r.returncode != 0 or not (out / "counts.f32").is_file():
            tail = "\n".join((r.stderr or r.stdout).strip().splitlines()[-5:])
            print(f"scnascope: {rds.name} not read: {tail}", file=sys.stderr)
            return None
        n, m = (int(v) for v in (out / "dims.txt").read_text().split())
        names = (out / "cells.txt").read_text().splitlines()
        if n != n_win:
            print(
                f"scnascope: {rds.name} has {n:,} windows, the table {n_win:,}; not used",
                file=sys.stderr,
            )
            return None
        C = np.fromfile(out / "counts.f32", dtype="<f4").reshape(m, n)
    rows = pd.Index(names).get_indexer(cols)
    X = np.full((len(cols), n), np.nan, dtype=np.float32)
    X[rows >= 0] = C[rows[rows >= 0]]
    return X


def _table(path: Path) -> Path | None:
    if path.is_file():
        return path if path.name.endswith(".tsv") else None
    for p in (path / TABLE, path / "epiAneufinder_results" / TABLE):
        if p.is_file():
            return p
    return None


def _detect(path: Path) -> bool:
    t = _table(path)
    if t is None:
        return False
    with open(t) as f:
        head = f.readline().split()
    return head[:3] == ["seq", "start", "end"] and any(h.startswith(PREFIX) for h in head[3:9])


def parse_bin_size(text: str | int | None) -> int | None:
    """``1Mb``, ``500kb``, ``2e6`` or a number of bp → bp (None: keep the windows)."""
    if text is None or text == "":
        return None
    if isinstance(text, int):
        return text
    t = str(text).strip().lower().replace(" ", "")
    for unit, k in (("mb", 1_000_000), ("kb", 1_000), ("bp", 1)):
        if t.endswith(unit):
            return int(round(float(t[: -len(unit)]) * k))
    return int(round(float(t)))


def merge_windows(
    bins: pd.DataFrame, state: np.ndarray, counts: np.ndarray | None, size: int
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray | None]:
    """Windows (sorted by chromosome and start) merged on a grid of `size` bp, per chromosome:
    counts summed, the majority call (ties are normal) per merged bin."""
    key = bins["chr"].astype(str) + ":" + ((bins["start"] - 1) // size).astype(str)
    codes, _ = pd.factorize(key)  # windows are sorted, so each merged bin is contiguous
    n = int(codes.max()) + 1 if len(codes) else 0
    g = bins.groupby(codes, sort=True)
    merged = pd.DataFrame(
        {
            "chr": g["chr"].first().to_numpy(),
            "start": g["start"].min().to_numpy(),
            "end": g["end"].max().to_numpy(),
            "n_windows": g.size().to_numpy(),
        }
    )
    votes = {}
    for v in (-1, 0, 1):
        m = np.zeros((state.shape[0], n), dtype=np.int32)
        np.add.at(m.T, codes, (state == v).T)
        votes[v] = m
    seen = votes[-1] + votes[0] + votes[1]
    best = np.maximum(votes[-1], np.maximum(votes[0], votes[1]))
    S = np.zeros_like(best, dtype=np.float32)
    S[(votes[1] == best) & (votes[-1] < best) & (votes[0] < best)] = 1
    S[(votes[-1] == best) & (votes[1] < best) & (votes[0] < best)] = -1
    S[seen == 0] = np.nan
    C = None
    if counts is not None:
        c = np.where(np.isfinite(counts), counts, 0).astype(np.float64)
        C = np.zeros((counts.shape[0], n))
        np.add.at(C.T, codes, c.T)
        C = C.astype(np.float32)
    return merged, S, C


def read_epianeufinder(
    path: Path,
    genome: str = "hg38",
    rscript: str | None = None,
    bin_size: str | int | None = None,
) -> CNVDataset:
    t = _table(path)
    if t is None:
        raise ValueError(f"no epiAneufinder {TABLE} in {path}")
    with open(t) as f:
        line = f.readline()
        head, first = line.split(), f.readline().split()
    sep = "\t" if "\t" in line else " "
    cols = [h for h in head[3:] if h.startswith(PREFIX)]
    # R's write.table adds row names: then data rows have one more field than the header
    names = ["_row", *head] if len(first) == len(head) + 1 else head
    raw = pd.read_csv(
        t,
        sep=sep,
        skiprows=1,
        header=None,
        names=names,
        dtype={"seq": str},
        na_values=["NA"],
        engine="c",
    )
    bins = pd.DataFrame(
        {
            "chr": raw["seq"].astype(str),
            "start": raw["start"].astype(np.int64),
            "end": raw["end"].astype(np.int64),
        }
    )
    order = np.lexsort((bins["start"].to_numpy(), bins["chr"].map(chr_sort_key).to_numpy()))
    bins = bins.iloc[order].reset_index(drop=True)
    calls = raw[cols].to_numpy(dtype=np.float32)[order].T  # cells x windows
    state = calls - 1  # 0/1/2 -> -1/0/+1
    fin = np.isfinite(state)
    n = np.maximum(fin.sum(axis=1), 1)
    cells = [c[len(PREFIX) :] for c in cols]
    obs = pd.DataFrame(
        {
            "cnv_burden": (np.where(fin, state != 0, False).sum(axis=1) / n),
            "gain_fraction": (np.where(fin, state > 0, False).sum(axis=1) / n),
            "loss_fraction": (np.where(fin, state < 0, False).sum(axis=1) / n),
        },
        index=cells,
    )
    counts = None
    if (t.parent / COUNTS).is_file():
        counts = _counts(t.parent / COUNTS, cols, len(bins), rscript)
        if counts is not None:
            counts = counts[:, order]
    # cells clustered with ward.D, as epiAneufinder's Karyogram, but on Manhattan distances
    # (as the scDNA tools here), on all of its windows (before any --bin-size merge)
    dendros = {}
    if 3 <= len(cells) <= 20000:
        dendros["total"] = ward_d_manhattan(state) | {
            "layers": ["state"],
            "method": "ward.D, Manhattan on epiAneufinder's calls (all its windows)",
        }
    size = parse_bin_size(bin_size)
    extra = {"dendrograms": dendros} if dendros else {}
    if size:
        window = int(np.median(bins["end"] - bins["start"] + 1))
        if size < 2 * window:
            raise ValueError(
                f"--bin-size {size:,} bp is not larger than epiAneufinder's windows "
                f"({window:,} bp); use at least {2 * window:,}"
            )
        bins, state, counts = merge_windows(bins, state, counts, size)
        extra |= {"bin_size": size, "window_size": window}
    layers = {"state": state.astype(np.float32)}
    sem = {"state": Semantics.TERNARY_STATE}
    if counts is not None:
        layers["read_counts"] = counts
        sem["read_counts"] = Semantics.COUNT
    root = t.parent.parent if t.parent.name == "epiAneufinder_results" else t.parent
    return CNVDataset.build(
        cells=obs,
        bins=bins,
        layers=layers,
        semantics=sem,
        modality=Modality.SCATAC,
        tool="epiAneufinder",
        genome=genome,
        default_layer="state",
        extra_uns={"source": root.name, **extra},
    )


#: a cell type needs this many cells to be its own reference; smaller ones use all cells
MIN_REFERENCE_CELLS = 20


def add_logr(ds: CNVDataset, by: str = "cell_type", min_cells: int = MIN_REFERENCE_CELLS) -> str:
    """Layer ``logr`` from ``read_counts``: per cell and bin, log2 of the cell's share of its
    fragments in the bin over the median share of the reference cells in that bin, centred on
    the cell's median bin. The reference is the cell's own cell type (``obs[by]``) when that
    has at least `min_cells` cells, otherwise all cells: chromatin accessibility differs between
    bins and between cell types, and dividing by the same bin in similar cells removes it.
    Negative counts (rare GC-correction artefacts) count as 0. Returns the reference used."""
    a = ds.adata
    if "read_counts" not in a.layers:
        return ""
    C = np.clip(np.nan_to_num(np.asarray(a.layers["read_counts"], dtype=np.float64)), 0, None)
    lib = C.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(lib > 0, (C + 0.5) / lib, np.nan)
    groups = np.zeros(a.n_obs, dtype=np.int64)
    names = ["all cells"]
    if by in a.obs:
        lab = a.obs[by].astype("object").to_numpy()
        valid = pd.Series(lab).notna().to_numpy()
        sizes = pd.Series(lab[valid]).value_counts()
        own = [k for k, n in sizes.items() if n >= min_cells]
        for k in own:
            names.append(str(k))
            groups[valid & (lab == k)] = len(names) - 1
    L = np.full(f.shape, np.nan)
    ref_all = np.nanmedian(f, axis=0)
    for gi in range(len(names)):
        rows = groups == gi
        if not rows.any():
            continue
        ref = ref_all if gi == 0 else np.nanmedian(f[rows], axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            L[rows] = np.log2(f[rows] / ref)
    with np.errstate(invalid="ignore"):
        L -= np.nanmedian(L, axis=1, keepdims=True)
    a.layers["logr"] = L.astype(np.float32)
    # MAPD: median |logR difference| of neighbouring bins within a chromosome (the noise)
    chrs = a.var["chr"].astype(str).to_numpy()
    same = chrs[1:] == chrs[:-1]
    d = np.where(same[None, :], np.abs(np.diff(L, axis=1)), np.nan)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        a.obs["mapd"] = np.nanmedian(d, axis=1)
    qc = dict(a.uns.get("qc", {}))
    if qc:
        qc["computed"] = sorted({*map(str, np.atleast_1d(qc.get("computed", []))), "mapd"})
        a.uns["qc"] = qc
    sem = dict(a.uns.get("semantics", {}))
    sem["logr"] = Semantics.CONTINUOUS_RATIO.value
    a.uns["semantics"] = sem
    own_types = names[1:]
    if own_types:
        ref = f"its own {by} ({len(own_types)} with >= {min_cells} cells; the others: all cells)"
    else:
        ref = "all cells"
    a.uns["logr_reference"] = ref
    return ref


def _read(path: Path, **opts) -> CNVDataset:
    return read_epianeufinder(path, **opts)


register(
    Reader(
        "epianeufinder",
        "scATAC",
        _detect,
        _read,
        "epiAneufinder output folder (epiAneufinder_results/results_table.tsv)",
    )
)
