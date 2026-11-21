"""Clustering of cells, computed from the canonical dataset.

Matches the ASCAT.sc Explorer and ASCAT.sc's R plots: Ward clustering (R's ``ward.D2``) of
the Manhattan distances between profiles, missing values set to -1. scipy's ``ward`` on a
precomputed distance vector applies the same Lance-Williams update as ``ward.D2``.
"""

from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from .model import CNVDataset, Semantics


def scipy_to_hclust(Z: np.ndarray) -> dict:
    """Convert a scipy linkage matrix to R's hclust ``merge``/``height``.

    hclust: row i merges two items; negative = leaf (-1-based leaf index), positive = the
    1-based row of an earlier merge. scipy: indices < n are leaves, n + i is merge row i.
    """
    n = Z.shape[0] + 1
    merge = np.empty((n - 1, 2), dtype=np.int64)
    for i in range(n - 1):
        for j in range(2):
            a = int(Z[i, j])
            merge[i, j] = -(a + 1) if a < n else (a - n + 1)
    return {"merge": merge, "height": Z[:, 2].astype(np.float64)}


def hclust_order(merge: np.ndarray) -> np.ndarray:
    """Leaf order of an hclust merge matrix (left child first), as row indices."""
    merge = np.asarray(merge, dtype=np.int64)
    if not len(merge):
        return np.zeros(1, dtype=np.int64)
    out: list[int] = []
    stack = [len(merge)]
    while stack:
        code = stack.pop()
        if code < 0:
            out.append(-code - 1)
            continue
        left, right = merge[code - 1]
        stack.append(int(right))
        stack.append(int(left))
    return np.asarray(out, dtype=np.int64)


def pool_columns(x: np.ndarray, max_cols: int) -> np.ndarray:
    """Means of blocks of neighbouring columns, so `x` has at most `max_cols` columns."""
    k = int(np.ceil(x.shape[1] / max_cols))
    if k <= 1:
        return x
    pad = (-x.shape[1]) % k
    xp = np.pad(x.astype(np.float64), ((0, 0), (0, pad)), constant_values=np.nan)
    with np.errstate(invalid="ignore"):
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmean(xp.reshape(x.shape[0], -1, k), axis=2)


def ward_manhattan(x: np.ndarray, max_cols: int | None = None) -> dict:
    """Ward.D2 on Manhattan distances of the rows of ``x`` (NaN → -1). With `max_cols`, a wide
    matrix is first averaged over blocks of neighbouring columns (thousands of cells on
    thousands of windows would otherwise take very long)."""
    if max_cols:
        x = pool_columns(x, max_cols)
    x = np.where(np.isnan(x), -1.0, x.astype(np.float64))
    Z = linkage(pdist(x, metric="cityblock"), method="ward")
    h = scipy_to_hclust(Z)
    h["order"] = leaves_list(Z).astype(np.int64)
    return h


#: wider matrices are averaged over blocks of neighbouring bins to this many columns before
#: clustering: 3,000 cells on 25,000 windows (epiAneufinder's 100 kb) took six minutes, 2,000
#: columns take seconds and keep every arm-level event
CLUSTER_COLS = 2000


def _manhattan_levels(x: np.ndarray) -> np.ndarray:
    """Manhattan distances between the rows of an integer-valued matrix (NaN → 0), by matrix
    products: |a - b| is the number of thresholds between consecutive levels that a and b fall
    on different sides of, and each threshold's disagreements are a Hamming distance."""
    x = np.nan_to_num(np.asarray(x, dtype=np.float64))
    levels = np.unique(x)
    d = np.zeros((len(x), len(x)))
    for lo, hi in zip(levels[:-1], levels[1:], strict=True):
        b = (x > lo).astype(np.float64)
        n = b.sum(axis=1)
        d += (hi - lo) * (n[:, None] + n[None, :] - 2.0 * (b @ b.T))
    np.maximum(d, 0, out=d)
    np.fill_diagonal(d, 0)
    return d


def ward_d_manhattan(x: np.ndarray) -> dict:
    """R's ``hclust(dist(x, "manhattan"), method = "ward.D")`` on the rows of ``x`` (NaN → 0),
    on every column. ward.D updates the distances themselves, which is scipy's ward
    (Lance-Williams on squared distances) given their square roots; the heights are squared
    back. A matrix with few integer levels (calls) uses exact matrix products, else pdist."""
    x = np.nan_to_num(np.asarray(x, dtype=np.float64))
    if len(np.unique(x)) <= 8 and np.allclose(x, np.round(x)):
        d = _manhattan_levels(x)[np.triu_indices(len(x), k=1)]
    else:
        d = pdist(x, metric="cityblock")
    Z = linkage(np.sqrt(d), method="ward")
    Z[:, 2] = Z[:, 2] ** 2
    h = scipy_to_hclust(Z)
    h["order"] = leaves_list(Z).astype(np.int64)
    return h


def cluster_cells(ds: CNVDataset, max_cells: int = 20000, force: bool = False) -> CNVDataset:
    """Add the ``total`` dendrogram (default layer) and, when the dataset has major and minor
    copy number, the ``allele`` dendrogram (major and minor stacked). Existing ones are kept
    unless ``force``. Skipped for fewer than 3 or more than ``max_cells`` cells.
    """
    n = ds.adata.n_obs
    if n < 3 or n > max_cells:
        return ds
    dendros = dict(ds.adata.uns.get("dendrograms", {}))
    if force or "total" not in dendros:
        layer = ds.default_layer
        x = ds.layer(layer)
        method = "ward.D2, Manhattan"
        if x.shape[1] > CLUSTER_COLS:
            method += f" (bins averaged to {CLUSTER_COLS} columns)"
        dendros["total"] = ward_manhattan(x, max_cols=CLUSTER_COLS) | {
            "layers": [layer],
            "method": method,
        }
    sem = ds.semantics
    if (
        (force or "allele" not in dendros)
        and sem.get("major_cn") is Semantics.ALLELE_CN
        and sem.get("minor_cn") is Semantics.ALLELE_CN
    ):
        x = np.hstack([ds.layer("major_cn"), ds.layer("minor_cn")])
        dendros["allele"] = ward_manhattan(x, max_cols=2 * CLUSTER_COLS) | {
            "layers": ["major_cn", "minor_cn"],
            "method": "ward.D2, Manhattan",
        }
    ds.adata.uns["dendrograms"] = dendros
    ds.validate()
    return ds
