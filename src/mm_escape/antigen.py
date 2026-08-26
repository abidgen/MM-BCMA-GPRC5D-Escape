"""Stage 08 — dual-antigen observed expression states, depth strata, nulls.

Every rule here is frozen in ``results/08_dual_antigen_escape/stage08_predeclaration.md``
and was written before any double-negative fraction existed.

The observed states are *measurements*, not biological antigen states: a zero is a zero,
and separating technical from biological zeros is the job of the depth-stratified null
and the sensitivity analyses, never of the call itself.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "ANTIGEN_GENES", "DN", "DP", "BCMA_ONLY", "GPRC5D_ONLY",
    "depth_ex_antigen", "observed_states", "quantile_edges", "assign_strata",
    "merge_sparse_strata", "stratified_expected_dn", "unconditioned_expected_dn",
    "permutation_null_dn", "hierarchical_bootstrap", "downsample_gene_counts",
]

#: Excluded from the depth metric so an antigen can never move a cell between strata.
ANTIGEN_GENES = ("TNFRSF17", "GPRC5D")

DP = "double_positive"
BCMA_ONLY = "BCMA_only"
GPRC5D_ONLY = "GPRC5D_only"
DN = "double_negative"


def depth_ex_antigen(total_counts, *antigen_counts):
    """Library size with the antigen genes' own UMIs removed.

    Using bare ``total_counts`` would let a perturbation of TNFRSF17 or GPRC5D shift a
    cell's depth stratum, which is exactly the circularity stage 07 was cleared of.
    """
    d = np.asarray(total_counts, dtype=np.float64).copy()
    for a in antigen_counts:
        d -= np.asarray(a, dtype=np.float64)
    return d


def observed_states(bcma_detected, gprc5d_detected):
    """Four observed expression states. Detection is binary and already applied."""
    b = np.asarray(bcma_detected, dtype=bool)
    g = np.asarray(gprc5d_detected, dtype=bool)
    out = np.empty(b.shape, dtype=object)
    out[b & g] = DP
    out[b & ~g] = BCMA_ONLY
    out[~b & g] = GPRC5D_ONLY
    out[~b & ~g] = DN
    return out


def quantile_edges(x, n_bins=5):
    """Quantile bin edges, de-duplicated.

    Heavy ties legitimately yield fewer than ``n_bins`` bins; that is reported rather
    than forced, because inventing an edge inside a tied mass would split identical
    cells into different technical regimes.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.size == 0:
        return np.array([], dtype=np.float64)
    qs = np.linspace(0.0, 1.0, n_bins + 1)
    return np.unique(np.quantile(x, qs))


def assign_strata(x, edges):
    """Left-closed / right-open bins; the top bin is closed. Returns 0-based indices."""
    x = np.asarray(x, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.float64)
    if edges.size < 2:
        return np.zeros(x.shape, dtype=np.int64)
    idx = np.searchsorted(edges, x, side="right") - 1
    return np.clip(idx, 0, edges.size - 2)


def merge_sparse_strata(strata, min_cells=20):
    """Merge strata below ``min_cells`` into the adjacent stratum by depth rank.

    The neighbour with more cells wins; ties go to the lower stratum. Iterates until
    every retained stratum clears the floor or a single stratum remains.
    """
    s = np.asarray(strata, dtype=np.int64).copy()
    while True:
        labels, counts = np.unique(s, return_counts=True)
        if labels.size <= 1 or counts.min() >= min_cells:
            return s
        j = int(np.argmin(counts))
        lab = labels[j]
        if j == 0:
            tgt = labels[1]
        elif j == labels.size - 1:
            tgt = labels[j - 1]
        else:
            lo, hi = labels[j - 1], labels[j + 1]
            tgt = hi if counts[j + 1] > counts[j - 1] else lo
        s[s == lab] = tgt


def stratified_expected_dn(bcma_neg, gprc5d_neg, strata):
    """Expected DN fraction under within-stratum independence.

    This is the closed form of the mean of the within-stratum permutation null, which is
    why the bootstrap can use it directly instead of permuting inside every replicate.
    It is a *technical baseline*, never a dropout correction: multiplying the marginals
    assumes exactly the independence the co-negativity test exists to interrogate.
    """
    b = np.asarray(bcma_neg, dtype=bool)
    g = np.asarray(gprc5d_neg, dtype=bool)
    s = np.asarray(strata)
    if b.size == 0:
        return np.nan
    total = 0.0
    for lab in np.unique(s):
        m = s == lab
        n = int(m.sum())
        total += n * b[m].mean() * g[m].mean()
    return total / b.size


def unconditioned_expected_dn(bcma_neg, gprc5d_neg):
    """Patient-wide P(BCMA=0) * P(GPRC5D=0), ignoring depth. Reported beside the
    stratified value; the gap between them is the depth artifact, quantified."""
    b = np.asarray(bcma_neg, dtype=bool)
    g = np.asarray(gprc5d_neg, dtype=bool)
    if b.size == 0:
        return np.nan
    return float(b.mean() * g.mean())


def permutation_null_dn(bcma_neg, gprc5d_neg, strata, n_perm=2000, seed=20260825):
    """Within-stratum permutation of the GPRC5D-zero label, BCMA held fixed.

    Preserves both marginals exactly within stratum and destroys only their
    within-stratum coupling. Returns the null DN fractions and a two-sided empirical p.
    """
    b = np.asarray(bcma_neg, dtype=bool)
    g = np.asarray(gprc5d_neg, dtype=bool)
    s = np.asarray(strata)
    rng = np.random.default_rng(seed)
    obs = float((b & g).mean()) if b.size else np.nan
    groups = [np.flatnonzero(s == lab) for lab in np.unique(s)]
    null = np.empty(n_perm, dtype=np.float64)
    for i in range(n_perm):
        gp = g.copy()
        for idx in groups:
            gp[idx] = rng.permutation(g[idx])
        null[i] = (b & gp).mean()
    p = (np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1) / (n_perm + 1)
    return null, float(p)


def hierarchical_bootstrap(sample_ids, values, n_boot=2000, seed=20260825, statistic=None):
    """Resample samples within a patient, then cells within each drawn sample.

    A flat cell bootstrap would treat sample-level batch variation as biological spread
    and report intervals that are too narrow. For a single-sample patient this
    degenerates to a cell bootstrap and cannot see sample-level variation at all — the
    caller is responsible for flagging that, since the arithmetic cannot.
    """
    sample_ids = np.asarray(sample_ids)
    values = np.asarray(values)
    rng = np.random.default_rng(seed)
    uniq = np.unique(sample_ids)
    by_sample = {u: np.flatnonzero(sample_ids == u) for u in uniq}
    if statistic is None:
        statistic = lambda v: float(np.mean(v))
    out = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        drawn = rng.choice(uniq, size=uniq.size, replace=True)
        picks = [rng.choice(by_sample[d], size=by_sample[d].size, replace=True) for d in drawn]
        out[i] = statistic(values[np.concatenate(picks)])
    return out


def downsample_gene_counts(total_umi, gene_counts, cap=10000, seed=20260825):
    """Exact subsample without replacement of one gene's counts to a library cap.

    ``Hypergeometric(N = total UMI, K = gene count, n = cap)``. Cells at or below the cap
    are returned unchanged — they are downsampled, never discarded, so the comparison
    stays over the same cells.
    """
    total = np.asarray(total_umi, dtype=np.int64)
    gene = np.asarray(gene_counts, dtype=np.int64)
    rng = np.random.default_rng(seed)
    out = gene.copy()
    hit = total > cap
    if not hit.any():
        return out
    ngood = gene[hit]
    nbad = total[hit] - ngood
    out[hit] = rng.hypergeometric(ngood=np.maximum(ngood, 0),
                                  nbad=np.maximum(nbad, 0), nsample=cap)
    return out
