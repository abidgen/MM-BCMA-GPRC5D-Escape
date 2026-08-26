"""Stage 09 — matched bulk RNA-seq (GSE223061) as orthogonal MARGINAL validation.

**Hard invariant, enforced by `tests/test_bulk.py`.** Bulk RNA can validate marginal
antigen abundance. It cannot determine whether the same individual cells are
simultaneously BCMA-negative and GPRC5D-negative: a tumour that is 50% BCMA+GPRC5D- plus
50% BCMA-GPRC5D+ shows healthy bulk expression of *both* genes while containing zero
dual-positive cells. Bulk destroys the joint distribution over cells by construction.

Nothing here may impute single-cell expression, convert a scRNA zero to a positive,
produce a corrected cell-level DN fraction, or scale the Stage-08 technical-zero floor.
Bulk is an independent measurement context, never a calibration equation.

The two bulk cohorts are different assays and are never pooled: MMRF is CD138+ sorted and
pairs with malignant-cell pseudobulk; WashU 1 is unsorted BMMC and pairs with whole-sample
pseudobulk.
"""
from __future__ import annotations

import gzip
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = [
    "ANTIGEN_ENSEMBL", "MATCHED_EXACT", "NOT_EVALUABLE",
    "antigen_transcript_map", "read_washu_bulk", "read_mmrf_bulk",
    "pseudobulk_cpm", "select_one_pair_per_patient", "spearman_or_not_evaluable",
]

ANTIGEN_ENSEMBL = {"TNFRSF17": "ENSG00000048462", "GPRC5D": "ENSG00000111291"}

MATCHED_EXACT = "MATCHED_EXACT"
NOT_EVALUABLE = "NOT_EVALUABLE"

#: A matched sample must clear the same >=20-cell floor stage 08 already uses for a
#: sample to carry its own estimate. Not a new constant.
MIN_CELLS_PER_MATCHED_SAMPLE = 20


def antigen_transcript_map(gtf_path, genes=tuple(ANTIGEN_ENSEMBL)):
    """Unversioned transcript ID -> gene symbol, for the antigens only.

    The WashU deposit is kallisto transcript-level output, so its TPM must be summed over
    a gene's transcripts before it is comparable to MMRF's gene-level TPM.
    """
    out = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if f[2] != "transcript":
                continue
            gn = re.search(r'gene_name "([^"]+)"', f[8])
            if gn and gn.group(1) in genes:
                tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
                out[tid.split(".")[0]] = gn.group(1)
    return out


def read_washu_bulk(tar_path, transcript_map):
    """Sum kallisto transcript TPM to gene for the antigens. Returns a dict, not cells."""
    with tarfile.open(tar_path, "r:gz") as tf:
        member = [x for x in tf.getnames() if x.endswith("abundance.tsv")][0]
        d = pd.read_csv(tf.extractfile(member), sep="\t")
    d["gene"] = d.target_id.str.split(".").str[0].map(transcript_map)
    agg = d.dropna(subset=["gene"]).groupby("gene").tpm.sum()
    return {"n_features": len(d), "n_runs": 1, "total_tpm": float(d.tpm.sum()),
            **{f"{g}_tpm": float(agg.get(g, np.nan)) for g in ANTIGEN_ENSEMBL}}


def read_mmrf_bulk(path):
    """Gene-level TPM. Multi-run files are AVERAGED, never summed.

    `MMRF_1686` stacks two sequencing runs in one file. TPM is already per-run normalised,
    so summing would double the apparent abundance.
    """
    try:
        d = pd.read_csv(path, sep="\t")
    except Exception:
        d = pd.DataFrame()
    if len(d) == 0:
        return {"n_features": 0, "n_runs": 0, "usable": False, "total_tpm": np.nan,
                **{f"{g}_tpm": np.nan for g in ANTIGEN_ENSEMBL}}
    runs = sorted(d.Sample_ID.unique())
    per = d[d.Gene_Name.isin(ANTIGEN_ENSEMBL)].groupby(["Sample_ID", "Gene_Name"]).TPM.sum().unstack()
    return {"n_features": int(len(d) / len(runs)), "n_runs": len(runs), "usable": True,
            "total_tpm": float(d.TPM.sum() / len(runs)),
            **{f"{g}_tpm": float(per[g].mean()) if g in per else np.nan for g in ANTIGEN_ENSEMBL}}


def pseudobulk_cpm(gene_counts, total_counts, mask):
    """Patient/sample-level abundance: 1e6 * sum(gene) / sum(total) over masked cells.

    The frozen Stage-09 primary metric. Aggregates sparsity rather than thresholding it,
    and needs no cutoff. Returns a scalar — never a per-cell value, so no bulk quantity
    can flow back onto a cell.
    """
    g = np.asarray(gene_counts)[mask]
    t = np.asarray(total_counts)[mask]
    if g.size == 0 or t.sum() == 0:
        return np.nan
    return float(1e6 * g.sum() / t.sum())


def select_one_pair_per_patient(matched, cell_counts, min_cells=MIN_CELLS_PER_MATCHED_SAMPLE):
    """One observation per biological patient: the earliest matched timepoint clearing
    the cell floor.

    Declared before any correlation was computed. Multiple samples from one patient are
    never counted as independent patients, and the timepoint is not chosen by which one
    agrees better with bulk.
    """
    df = matched.copy()
    df["n_cells"] = df.sc_sample.map(cell_counts).fillna(0).astype(int)
    df = df[df.n_cells >= min_cells]
    return (df.sort_values(["sc_patient", "timepoint"])
              .groupby("sc_patient", as_index=False).first())


def spearman_or_not_evaluable(x, y, n_min=5):
    """Spearman rho, or an explicit NOT_EVALUABLE at small n.

    An underpowered subgroup is reported as not evaluable, never as "no relationship" —
    the p-value is descriptive support only and is not interpreted at small matched n.
    """
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < n_min:
        return {"n": int(ok.sum()), "rho": np.nan, "p": np.nan,
                "status": f"{NOT_EVALUABLE} (n<{n_min})"}
    rho, p = spearmanr(x[ok], y[ok])
    return {"n": int(ok.sum()), "rho": float(rho), "p": float(p), "status": "evaluable"}
