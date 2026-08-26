"""Stage 10 — is the DN population a coherent transcriptional STATE?

Design frozen in ``results/10_dn_coherence/stage10_design.md``.

**The vocabulary rule is load-bearing.** Transcriptional coherence licenses
"escape-associated transcriptional state" and nothing stronger. "Subclone" requires
independent CNV/genetic support, which this project does not have — Stage 07's CNV
inference failed its donor negative control and is frozen NOT_EVALUABLE, and that failed
method is not silently reused here.

**The antigen-circularity rule is enforced structurally, not by convention.** The antigen
columns are dropped from the matrix *before normalization*, so nothing downstream — the
library-size factor included — can carry information about them. Selecting DN cells by
antigen status is legitimate; letting antigen expression into the feature matrix that then
"discovers" DN structure is not.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

__all__ = [
    "ANTIGEN_FEATURES", "MIN_GROUP_CELLS", "MIN_PATIENT_CELLS", "N_DEPTH_BINS",
    "MAX_DEPTH_BINS", "adaptive_depth_bins",
    "SUPPORTED", "NOT_SUPPORTED", "NOT_EVALUABLE",
    "CNV_NOT_EVALUABLE", "drop_antigen_features", "depth_bins",
    "knn_dn_fraction", "morans_i", "best_cluster_enrichment",
    "depth_stratified_permutation", "depth_matched_indices", "coherence_state",
    "STRUCTURE_SUPPORTED", "STRUCTURE_NOT_SUPPORTED", "STRUCTURE_NOT_EVALUABLE",
    "STATE_SUPPORTED", "STATE_NOT_SUPPORTED", "STATE_NOT_EVALUABLE",
    "CNV_SUBCLONE_SUPPORTED", "level1_structure_state", "level2_state",
    "licensed_language",
]

#: Never features. Labels only.
ANTIGEN_FEATURES = ("TNFRSF17", "GPRC5D")

MIN_GROUP_CELLS = 20      # frozen stage-08 floor for a unit to carry its own estimate
MIN_PATIENT_CELLS = 100   # frozen stage-08 minimum-n rule
N_DEPTH_BINS = 5
MAX_DEPTH_BINS = 10

SUPPORTED = "DN_COHERENCE_SUPPORTED"
NOT_SUPPORTED = "DN_COHERENCE_NOT_SUPPORTED"
NOT_EVALUABLE = "DN_COHERENCE_NOT_EVALUABLE"

#: The only CNV state this project may emit. Not "not supported" — the assay failed its
#: negative control, so it is unevaluable, and an underpowered null is not a negative.
CNV_NOT_EVALUABLE = "CNV_SUBCLONE_NOT_EVALUABLE"

#: Emitted ONLY from independently validated CNV evidence, which this project does not
#: have. There is deliberately no `CNV_SUBCLONE_NOT_SUPPORTED`: the assay failed its
#: donor negative control, so it is unevaluable, and an underpowered null is not a negative.
CNV_SUBCLONE_SUPPORTED = "CNV_SUBCLONE_SUPPORTED"

# Level 1 — spatial/transcriptional-space enrichment. Licenses "non-random DN
# organization" and NOTHING more; it is not an escape-associated state.
STRUCTURE_SUPPORTED = "DN_STRUCTURE_SUPPORTED"
STRUCTURE_NOT_SUPPORTED = "DN_STRUCTURE_NOT_SUPPORTED"
STRUCTURE_NOT_EVALUABLE = "DN_STRUCTURE_NOT_EVALUABLE"

# Level 2 — reproducible non-antigen program. Licenses "escape-associated
# transcriptional state".
STATE_SUPPORTED = "DN_STATE_SUPPORTED"
STATE_NOT_SUPPORTED = "DN_STATE_NOT_SUPPORTED"
STATE_NOT_EVALUABLE = "DN_STATE_NOT_EVALUABLE"


def drop_antigen_features(matrix, gene_names, antigens=ANTIGEN_FEATURES):
    """Remove the antigen columns *before* any normalization.

    Order matters and is the whole point: normalising first would fold antigen counts into
    every cell's size factor, so perturbing an antigen would shift every other gene's
    normalised value and the "independent" coherence result would not be independent.
    """
    keep = [i for i, g in enumerate(gene_names) if g not in set(antigens)]
    return matrix[:, keep], [gene_names[i] for i in keep]


def depth_bins(depth, n_bins=N_DEPTH_BINS):
    """Within-patient quantile bins of an antigen-independent depth metric."""
    depth = np.asarray(depth, dtype=float)
    if depth.size == 0:
        return np.zeros(0, dtype=np.int64)
    edges = np.unique(np.quantile(depth, np.linspace(0, 1, n_bins + 1)))
    if edges.size < 2:
        return np.zeros(depth.size, dtype=np.int64)
    return np.clip(np.searchsorted(edges, depth, side="right") - 1, 0, edges.size - 2)


def adaptive_depth_bins(depth, min_per_bin=MIN_GROUP_CELLS, max_bins=MAX_DEPTH_BINS):
    """As many quantile depth bins as keep `min_per_bin` cells each, capped at `max_bins`.

    Amended before any patient data was analysed, because a synthetic diagnostic showed
    five quintiles leave enough residual within-bin depth variation for a purely
    depth-driven label to still read as coherent (p = 0.044); ten bins gave 0.19. The rule
    is tied to the existing frozen 20-cell floor rather than to that number, and it costs
    power against real signal — the conservative direction.
    """
    depth = np.asarray(depth, dtype=float)
    n_bins = int(np.clip(depth.size // max(min_per_bin, 1), 2, max_bins))
    return depth_bins(depth, n_bins)


def knn_dn_fraction(labels, neighbors):
    """Mean fraction of a DN cell's neighbours that are also DN.

    `neighbors` is an (n_cells, k) integer array of neighbour indices.
    """
    labels = np.asarray(labels, dtype=bool)
    if labels.sum() == 0:
        return np.nan
    return float(labels[neighbors[labels]].mean())


def morans_i(labels, graph):
    """Moran's I of the DN label over a (sparse, symmetric) neighbour graph."""
    x = np.asarray(labels, dtype=float)
    n = x.size
    w = sp.csr_matrix(graph)
    z = x - x.mean()
    denom = (z ** 2).sum()
    s0 = w.sum()
    if denom == 0 or s0 == 0:
        return np.nan
    return float((n / s0) * (z @ (w @ z)) / denom)


def best_cluster_enrichment(labels, clusters):
    """Largest DN over-representation in any local cluster, as a ratio to the patient rate.

    Reported, but never sufficient on its own: one local cluster enriched for DN cells is
    not automatically a biological state.
    """
    labels = np.asarray(labels, dtype=bool)
    clusters = np.asarray(clusters)
    base = labels.mean()
    if base in (0.0, 1.0):
        return np.nan
    best = 0.0
    for c in np.unique(clusters):
        m = clusters == c
        if m.sum() < MIN_GROUP_CELLS:
            continue
        best = max(best, labels[m].mean() / base)
    return float(best)


def depth_stratified_permutation(labels, bins, statistic, n_perm=1000, seed=20260825):
    """Shuffle DN labels **within depth strata** and recompute the statistic.

    Shallow cells are both more likely to be DN and more likely to sit together in
    low-dimensional space, so an unconditioned permutation reports "coherence" from library
    size alone — an artifact pointing in exactly the direction this project hopes to find.
    Stratifying preserves the depth-label coupling and destroys only the
    transcription-label coupling.

    Returns (null distribution, one-sided empirical p for observed > null).
    """
    labels = np.asarray(labels, dtype=bool)
    bins = np.asarray(bins)
    rng = np.random.default_rng(seed)
    observed = statistic(labels)
    groups = [np.flatnonzero(bins == b) for b in np.unique(bins)]
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm = labels.copy()
        for idx in groups:
            perm[idx] = rng.permutation(labels[idx])
        null[i] = statistic(perm)
    ok = np.isfinite(null)
    p = ((np.sum(null[ok] >= observed) + 1) / (ok.sum() + 1)) if ok.any() else np.nan
    return null, float(p)


def depth_matched_indices(is_dn, bins, seed=20260825):
    """Per depth bin, take min(n_dn, n_pos) from each group.

    Frozen before results as the DE handling, so a depth difference between DN and
    antigen-positive cells cannot masquerade as a transcriptional program.
    """
    is_dn = np.asarray(is_dn, dtype=bool)
    bins = np.asarray(bins)
    rng = np.random.default_rng(seed)
    dn_out, pos_out = [], []
    for b in np.unique(bins):
        d = np.flatnonzero((bins == b) & is_dn)
        p = np.flatnonzero((bins == b) & ~is_dn)
        k = min(d.size, p.size)
        if k == 0:
            continue
        dn_out.append(rng.choice(d, k, replace=False))
        pos_out.append(rng.choice(p, k, replace=False))
    if not dn_out:
        return np.array([], dtype=int), np.array([], dtype=int)
    return np.sort(np.concatenate(dn_out)), np.sort(np.concatenate(pos_out))


def coherence_state(primary, sensitivity, repeated_sample_status=None, alpha=0.05):
    """Collapse the two denominators into one three-state call.

    `primary`/`sensitivity` are dicts with `evaluable` and `perm_p` (the depth-stratified
    p). Requiring both denominators is the same discipline Stage 08 fixed: if coherence
    appears under only one, that is instability, not a result to be selected.
    """
    if not primary.get("evaluable") or not sensitivity.get("evaluable"):
        return NOT_EVALUABLE
    if repeated_sample_status == "discordant":
        return NOT_SUPPORTED
    pp, ps = primary.get("perm_p", np.nan), sensitivity.get("perm_p", np.nan)
    if np.isfinite(pp) and np.isfinite(ps) and pp < alpha and ps < alpha:
        return SUPPORTED
    return NOT_SUPPORTED


def level1_structure_state(coherence_state_value):
    """Map the original Level-1 call onto explicit structure vocabulary.

    The original label was `DN_COHERENCE_SUPPORTED`, which reads as though it licensed an
    escape-associated *state*. It never did: it rested on non-random local organization
    only. The numbers are unchanged; only the word is corrected.
    """
    return {SUPPORTED: STRUCTURE_SUPPORTED,
            NOT_SUPPORTED: STRUCTURE_NOT_SUPPORTED,
            NOT_EVALUABLE: STRUCTURE_NOT_EVALUABLE}[coherence_state_value]


def level2_state(evaluable, reproducible_program_hits, repeated_sample_status=None):
    """Level 2 from PROGRAM evidence only.

    `reproducible_program_hits` are programs that were cohort-reproducible under both
    denominators *and* run the same direction in this patient. Neighbourhood enrichment
    and Moran's I are Level-1 quantities and cannot reach this function at all — which is
    the point: Level 2 may never be emitted from spatial organization.
    """
    if not evaluable:
        return STATE_NOT_EVALUABLE
    if repeated_sample_status == "discordant":
        return STATE_NOT_SUPPORTED
    return STATE_SUPPORTED if list(reproducible_program_hits) else STATE_NOT_SUPPORTED


def licensed_language(level1, level2, level3=CNV_NOT_EVALUABLE):
    """What each evidence level permits you to write. Nothing stronger is allowed."""
    if level3 == CNV_SUBCLONE_SUPPORTED:
        return "escape-associated subclone"
    if level2 == STATE_SUPPORTED:
        return "escape-associated transcriptional state"
    if level1 == STRUCTURE_SUPPORTED:
        return "non-random DN organization"
    return "none"
