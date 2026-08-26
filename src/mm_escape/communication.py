"""Stage 11 — exploratory immune context. Design frozen in
``results/11_immune_context/stage11_design.md``.

**Stage 11 asks whether immune context accompanies the frozen DN measurement and coherence
phenotypes. It does not decide whether those phenotypes are real.** Nothing here may
rescue, upgrade or downgrade a patient, modify a frozen tier, or emit a scalar risk score.

Stages 08 and 10 both produced strong apparent biology out of sequencing depth, so the
default assumption here is that an unadjusted immune association is a technical one until
cohort, depth and cell-count adjustment says otherwise. The unadjusted estimate is always
reported beside the adjusted one, because the gap between them is the informative part.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "IMMUNE_CATEGORIES", "MIN_IMMUNE_CELLS", "NOT_EVALUABLE", "CONFOUNDERS",
    "MEASUREMENT_PREDICTORS", "LR_CANDIDATES",
    "LR_SENDERS", "MIN_SENDER_CELLS",
    "clr", "ols_association", "benjamini_hochberg", "within_cohort_spearman",
    "stream_gene_counts", "pseudobulk_cpm",
]

#: Frozen from the accepted stage-06 endpoint plus the cluster-23 TRBC-context revision.
#: `cytotoxic_mixed` is never folded into T or NK.
IMMUNE_CATEGORIES = ("Tcell", "NK_core", "cytotoxic_mixed", "Myeloid", "Bcell",
                     "HSPC", "Erythroid")

#: The frozen stage-08 MIN_PATIENT_CELLS, reused rather than re-derived.
MIN_IMMUNE_CELLS = 100

NOT_EVALUABLE = "NOT_EVALUABLE"

#: Fixed before any association was computed. Cohort is required, not optional.
CONFOUNDERS = ("cohort", "log_depth", "log_n_immune", "log_n_samples")

#: All four are reported. None is selected for correlating best.
MEASUREMENT_PREDICTORS = ("obs_dn_primary", "obs_dn_sensitivity", "excess_dn",
                          "enr_cohortbins")

#: Predeclared candidate ligand->receptor pairs, sender immune -> receiver clone plasma.
#: Discovery-only context; a full-interactome screen is deliberately not run at n ~ 32.
LR_CANDIDATES = (
    ("PRF1", "None"), ("GZMB", "None"), ("GZMA", "F2R"), ("GZMK", "F2R"),
    ("IFNG", "IFNGR1"), ("IFNG", "IFNGR2"),
    ("FASLG", "FAS"), ("TNFSF10", "TNFRSF10A"), ("TNFSF10", "TNFRSF10B"),
    ("PDCD1", "CD274"), ("CTLA4", "CD86"), ("TIGIT", "PVR"),
    ("KLRK1", "MICA"), ("KLRK1", "MICB"), ("KLRK1", "ULBP2"),
    ("TNFSF13", "TNFRSF13B"), ("TNFSF13B", "TNFRSF13B"),
    ("IL6", "IL6R"), ("CXCL12", "CXCR4"),
)


def clr(fractions, pseudocount=1e-6):
    """Centred log-ratio of a compositional row.

    Immune fractions sum to one, so one lineage rising forces the others down. Treating
    each as an independent quantity is what makes naive per-type tests anticonservative;
    the CLR removes the constraint before any model sees the data.
    """
    x = np.asarray(fractions, dtype=float) + pseudocount
    log_x = np.log(x)
    return log_x - log_x.mean(axis=-1, keepdims=True)


def _design(predictor, confounders):
    cols = [np.ones(len(predictor)), np.asarray(predictor, dtype=float)]
    for c in confounders:
        c = np.asarray(c)
        if c.dtype.kind in "OUS":                      # categorical -> dummies, drop first
            levels = sorted(set(c.tolist()))[1:]
            cols.extend([(c == lv).astype(float) for lv in levels])
        else:
            cols.append(c.astype(float))
    return np.column_stack(cols)


def ols_association(feature, predictor, confounders=(), alpha=0.05):
    """OLS coefficient on `predictor`, adjusted for `confounders`. Patient is the unit.

    Returns coefficient, standard error, t, p, confidence interval, n and residual df.
    Deliberately plain: at n ~ 32 the effect size and its interval carry the message, and
    a p-value is supporting detail rather than a finding.
    """
    from scipy import stats

    y = np.asarray(feature, dtype=float)
    ok = np.isfinite(y) & np.isfinite(np.asarray(predictor, dtype=float))
    for c in confounders:
        c = np.asarray(c)
        if c.dtype.kind not in "OUS":
            ok &= np.isfinite(c.astype(float))
    y = y[ok]
    X = _design(np.asarray(predictor, dtype=float)[ok],
                [np.asarray(c)[ok] for c in confounders])
    n, k = X.shape
    df = n - k
    if df <= 0 or n < 5:
        return {"n": int(n), "coef": np.nan, "se": np.nan, "t": np.nan, "p": np.nan,
                "ci_lo": np.nan, "ci_hi": np.nan, "df": int(df), "status": NOT_EVALUABLE}
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    sigma2 = (resid @ resid) / df
    xtx_inv = np.linalg.pinv(X.T @ X)
    se = float(np.sqrt(sigma2 * xtx_inv[1, 1]))
    coef = float(beta[1])
    t = coef / se if se > 0 else np.nan
    p = float(2 * stats.t.sf(abs(t), df)) if np.isfinite(t) else np.nan
    crit = stats.t.ppf(1 - alpha / 2, df)
    return {"n": int(n), "coef": coef, "se": se, "t": float(t), "p": p,
            "ci_lo": coef - crit * se, "ci_hi": coef + crit * se,
            "df": int(df), "status": "evaluable"}


def benjamini_hochberg(pvals):
    """BH-adjusted p-values, NaNs passed through."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if ok.size == 0:
        return out
    order = ok[np.argsort(p[ok])]
    m = order.size
    adj = p[order] * m / np.arange(1, m + 1)
    out[order] = np.minimum.accumulate(adj[::-1])[::-1].clip(0, 1)
    return out


def within_cohort_spearman(feature, predictor, cohort, n_min=5):
    """Per-cohort Spearman, or an explicit NOT_EVALUABLE.

    An underpowered subgroup is reported as not evaluable, never as "no relationship".
    """
    from scipy.stats import spearmanr

    feature = np.asarray(feature, dtype=float)
    predictor = np.asarray(predictor, dtype=float)
    cohort = np.asarray(cohort)
    out = {}
    for c in sorted(set(cohort.tolist())):
        m = (cohort == c) & np.isfinite(feature) & np.isfinite(predictor)
        if m.sum() < n_min:
            out[c] = {"n": int(m.sum()), "rho": np.nan, "status": f"{NOT_EVALUABLE} (n<{n_min})"}
            continue
        rho, p = spearmanr(feature[m], predictor[m])
        out[c] = {"n": int(m.sum()), "rho": float(rho), "p": float(p), "status": "evaluable"}
    return out


#: Sender categories for the communication arm. A category must reach `MIN_SENDER_CELLS`
#: in a patient before it is allowed to carry that patient's ligand estimate.
LR_SENDERS = ("Tcell", "NK_core", "cytotoxic_mixed", "Myeloid")

#: The frozen stage-08 `MIN_GROUP_CELLS`, reused rather than re-derived.
MIN_SENDER_CELLS = 20


def stream_gene_counts(h5ad_path, genes, layer="counts", chunk=20000):
    """Per-cell counts for `genes`, plus each cell's total over the stored gene space.

    Reads the CSR layer in row blocks rather than materialising it. The stage-05 object is
    172,940 x 32,991 and loading it costs ~20 GB; this stage needs about twenty columns.

    The returned `total_counts` is the row sum **over the intersected gene space**, which is
    what stage 08 normalised against. It is deliberately not `obs["total_counts"]`, which
    was computed at QC time over each sample's full Cell Ranger reference and therefore runs
    a few counts higher — enough to move a per-patient median.
    """
    import h5py
    import pandas as pd
    from anndata.io import read_elem

    genes = list(genes)
    with h5py.File(h5ad_path, "r") as f:
        var_index = read_elem(f["var"]).index.astype(str)
        obs_index = read_elem(f["obs"]).index.astype(str)
        missing = [g for g in genes if g not in set(var_index)]
        if missing:
            raise KeyError(f"genes absent from {h5ad_path}: {missing}")
        pos = {g: int(np.flatnonzero(var_index == g)[0]) for g in genes}
        col_of = {c: j for j, (g, c) in enumerate(pos.items())}

        grp = f["layers"][layer] if layer else f["X"]
        indptr = grp["indptr"][:]
        data, indices = grp["data"], grp["indices"]
        n = len(indptr) - 1
        out = np.zeros((n, len(genes)))
        total = np.zeros(n)
        wanted = np.array(sorted(col_of), dtype=indices.dtype)
        for start in range(0, n, chunk):
            stop = min(start + chunk, n)
            lo, hi = indptr[start], indptr[stop]
            d, ix = data[lo:hi], indices[lo:hi]
            rel = indptr[start:stop + 1] - lo
            rows = np.repeat(np.arange(start, stop), np.diff(rel))
            np.add.at(total, rows, d)
            m = np.isin(ix, wanted)
            if m.any():
                np.add.at(out, (rows[m], np.array([col_of[c] for c in ix[m]])), d[m])

    frame = pd.DataFrame(out, columns=genes, index=pd.Index(obs_index, name="cell_id"))
    frame["total_counts"] = total
    return frame


def pseudobulk_cpm(counts, genes, total_col="total_counts"):
    """Counts per million over a pooled group of cells.

    Pooled, not a mean of per-cell rates: a mean over cells weights a 300-UMI cell and a
    20,000-UMI cell equally, which in this cohort means weighting a WashU cell like an MMRF
    one. Returns NaN for an empty group rather than zero — absent is not the same as zero.
    """
    import pandas as pd

    total = float(counts[total_col].sum())
    if len(counts) == 0 or total <= 0:
        return pd.Series(np.nan, index=list(genes))
    return counts[list(genes)].sum() / total * 1e6
