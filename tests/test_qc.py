"""Stage 04 QC — thresholds, outlier calling, doublets, checkpointing.

Two tiers, per `tests/conftest.py`: the threshold arithmetic and the flag logic are
exercised on synthetic frames with no deposit at all; anything that touches a real
matrix, and the scDblFinder bridge, sit behind `requires_data`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from mm_escape import qc

from conftest import requires_data, requires_r


# ---------------------------------------------------------------------------
# The MAD rule — no data needed
# ---------------------------------------------------------------------------

def test_mad_outlier_flags_both_tails():
    rng = np.random.default_rng(0)
    values = np.concatenate([rng.normal(10.0, 1.0, 200), [-20.0, 40.0]])
    flagged = qc.mad_outlier(values, n_mads=5.0)
    assert flagged[-2:].all()
    assert not flagged[:200].any()


def test_mad_outlier_uses_the_raw_unscaled_mad():
    """Raw MAD, not 1.4826*MAD. Mixing the two moves every threshold by a third."""
    values = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])        # median 0, raw MAD 1
    # At 1.5 MADs the raw band is |v| > 1.5, which flags the two ends. Under the
    # scaled MAD (1.4826) the same call would be |v| > 2.22 and would flag nothing,
    # so this asserts the convention rather than merely exercising the function.
    assert list(qc.mad_outlier(values, n_mads=1.5)) == [
        True, False, False, False, True
    ]
    assert not qc.mad_outlier(values, n_mads=1.5 * 1.4826).any()


def test_mad_outlier_on_a_constant_metric_flags_nothing():
    """A zero MAD makes every non-median point infinitely deviant. Filter nothing."""
    assert not qc.mad_outlier(np.full(50, 7.0), n_mads=5.0).any()
    # One deviant point among constants still gives MAD 0 — and must not wipe the set.
    values = np.append(np.full(50, 7.0), 99.0)
    assert not qc.mad_outlier(values, n_mads=5.0).any()


def test_mad_thresholds_reports_the_interval_and_both_tails():
    frame = pd.DataFrame({"m": [0.0, 1.0, 2.0, 3.0, 4.0, 50.0]})
    table = qc.mad_thresholds(frame, ["m"], n_mads=2.0)
    row = table.iloc[0]
    assert row["metric"] == "m"
    assert row["median"] == 2.5
    assert row["mad"] == 1.5
    assert row["lower"] == pytest.approx(-0.5)
    assert row["upper"] == pytest.approx(5.5)
    assert row["n_below"] == 0
    assert row["n_above"] == 1


def test_default_mad_counts_are_the_standard_five():
    """Not tightened ahead of the data — see the constant's own note."""
    assert set(qc.DEFAULT_N_MADS.values()) == {5.0}
    assert set(qc.DEFAULT_N_MADS) == set(qc.MAD_METRICS)


# ---------------------------------------------------------------------------
# flag_outliers — synthetic AnnData, still no deposit
# ---------------------------------------------------------------------------

def _synthetic(n_per_group=400, seed=0):
    """Two cohorts differing ~2x in depth, like MMRF vs WashU 1."""
    rng = np.random.default_rng(seed)
    deep = rng.normal(9.0, 0.3, n_per_group)
    shallow = rng.normal(8.0, 0.3, n_per_group)
    obs = pd.DataFrame({
        "cohort": ["MMRF"] * n_per_group + ["WU1"] * n_per_group,
        "log1p_total_counts": np.concatenate([deep, shallow]),
        "log1p_n_genes_by_counts": np.concatenate([deep - 1, shallow - 1]),
        "pct_counts_in_top_20_genes": rng.normal(30, 3, 2 * n_per_group),
        "pct_counts_mt": rng.normal(8, 2, 2 * n_per_group),
        "n_genes_by_counts": np.exp(np.concatenate([deep - 1, shallow - 1])),
    }, index=[f"c{i}" for i in range(2 * n_per_group)])
    return AnnData(
        X=np.zeros((len(obs), 0), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index([], name="g")),
    )


def test_flag_outliers_requires_metrics_first():
    adata = AnnData(X=np.zeros((3, 0), dtype=np.float32))
    with pytest.raises(ValueError, match="add_qc_metrics"):
        qc.flag_outliers(adata)


def test_flag_outliers_writes_every_flag_but_unions_only_the_filters():
    adata = _synthetic()
    qc.flag_outliers(adata)
    for column in [*qc.ALL_FLAGS, "outlier"]:
        assert adata.obs[column].dtype == bool
    union = np.logical_or.reduce(
        [adata.obs[c].to_numpy() for c in qc.DEFAULT_FILTERS]
    )
    assert (adata.obs["outlier"].to_numpy() == union).all()
    assert adata.uns["qc_filters"] == list(qc.DEFAULT_FILTERS)


def test_top20_is_computed_but_does_not_delete_cells():
    """The stage's one real departure from the generic recipe.

    In myeloma marrow a library dominated by a few transcripts is a plasma cell's
    normal state, so a MAD filter on `pct_counts_in_top_20_genes` deletes
    antigen-POSITIVE malignant cells and inflates the escape fraction. The flag is
    still computed — it is an ambient-Ig handle stage 08 needs — but it is not in
    `DEFAULT_FILTERS`.
    """
    assert "outlier_top20" in qc.ALL_FLAGS
    assert "outlier_top20" not in qc.DEFAULT_FILTERS

    adata = _synthetic()
    adata.obs.loc["c0", "pct_counts_in_top_20_genes"] = 95.0
    qc.flag_outliers(adata, group_key=None)
    assert adata.obs.loc["c0", "outlier_top20"], "the flag must still be computed"
    assert not adata.obs.loc["c0", "outlier"], "but it must not delete the cell"

    # Opting back in is possible and explicit, for a sensitivity re-run.
    qc.flag_outliers(adata, group_key=None, filters=(*qc.DEFAULT_FILTERS, "outlier_top20"))
    assert adata.obs.loc["c0", "outlier"]


def test_flag_outliers_rejects_an_unknown_filter():
    with pytest.raises(ValueError, match="unknown filter"):
        qc.flag_outliers(_synthetic(), filters=("outlier_counts", "outlier_nonsense"))


def test_per_cohort_thresholds_differ_from_pooled_ones():
    """The stage's central decision: a pooled MAD punishes the shallow cohort.

    With two cohorts a full log1p unit apart, pooling widens the median and cuts
    into the shallow cohort's own distribution. Per-cohort thresholds must sit
    inside each cohort and must not be equal to each other.
    """
    adata = _synthetic()
    per_cohort = qc.flag_outliers(adata.copy(), group_key="cohort")
    counts = per_cohort.query("metric == 'log1p_total_counts'").set_index("group")
    assert set(counts.index) == {"MMRF", "WU1"}
    assert counts.loc["MMRF", "median"] > counts.loc["WU1", "median"] + 0.5
    assert counts.loc["MMRF", "lower"] > counts.loc["WU1", "lower"]

    pooled = qc.flag_outliers(adata.copy(), group_key=None)
    pooled_row = pooled.query("metric == 'log1p_total_counts'").iloc[0]
    # The pooled median sits between the two cohorts, so it describes neither.
    assert counts.loc["WU1", "median"] < pooled_row["median"] < counts.loc["MMRF", "median"]
    # And the pooled MAD is inflated by the between-cohort gap.
    assert pooled_row["mad"] > counts["mad"].max()


def test_mitochondrial_rule_is_one_sided():
    """A cell with unusually FEW mitochondrial reads is not low quality."""
    adata = _synthetic()
    adata.obs.loc["c0", "pct_counts_mt"] = 0.0
    adata.obs.loc["c1", "pct_counts_mt"] = 60.0
    table = qc.flag_outliers(adata, group_key=None)
    assert not adata.obs.loc["c0", "outlier_mt"]
    assert adata.obs.loc["c1", "outlier_mt"]
    row = table.query("metric == 'pct_counts_mt'").iloc[0]
    assert row["lower"] == -np.inf
    assert row["n_below"] == 0


def test_mitochondrial_cap_binds_and_says_so():
    adata = _synthetic()
    table = qc.flag_outliers(adata, group_key=None, mt_max_pct=9.0)
    row = table.query("metric == 'pct_counts_mt'").iloc[0]
    assert row["upper"] == 9.0
    assert row["mt_bound_by"] == "cap"
    loose = qc.flag_outliers(_synthetic(), group_key=None, mt_max_pct=100.0)
    assert loose.query("metric == 'pct_counts_mt'").iloc[0]["mt_bound_by"] == "mad"


def test_min_genes_is_an_absolute_floor_on_top_of_the_mad_rule():
    adata = _synthetic()
    adata.obs.loc["c0", "n_genes_by_counts"] = 10.0
    qc.flag_outliers(adata, group_key=None, min_genes=200)
    assert adata.obs.loc["c0", "outlier_min_genes"]
    assert adata.obs.loc["c0", "outlier"]


def test_flag_outliers_rejects_a_missing_group_key():
    adata = _synthetic()
    del adata.obs["cohort"]
    with pytest.raises(ValueError, match="cohort"):
        qc.flag_outliers(adata, group_key="cohort")


def test_cohort_thresholds_matches_flag_outliers_on_the_same_obs():
    """`cohort_thresholds` exists so the stage never concatenates matrices; it must
    be the same rule, not a second implementation of it."""
    adata = _synthetic()
    reference = adata.copy()
    direct = qc.flag_outliers(reference, group_key="cohort")
    flags, table = qc.cohort_thresholds(adata.obs.copy(), group_key="cohort")

    pd.testing.assert_frame_equal(direct, table)
    assert len(flags) == adata.n_obs
    assert list(flags.columns[-1:]) == ["outlier"]
    for column in flags.columns:
        assert (flags[column].to_numpy() == reference.obs[column].to_numpy()).all()


# ---------------------------------------------------------------------------
# Against the real deposit
# ---------------------------------------------------------------------------

@requires_data
def test_qc_metrics_on_a_real_sample(samples):
    adata = samples["BM4"].copy()
    qc.add_qc_metrics(adata)
    for column in ("pct_counts_in_top_20_genes", "pct_counts_mt", "pct_counts_ribo",
                   "pct_counts_hb", "log1p_total_counts",
                   "log1p_n_genes_by_counts"):
        assert column in adata.obs
    # 13 protein-coding mitochondrial genes in GRCh38, in both reference builds.
    assert adata.var["mt"].sum() == 13
    assert 0 <= adata.obs["pct_counts_mt"].max() <= 100


@requires_data
def test_haemoglobin_match_does_not_catch_HBEGF(samples):
    """The trailing character class matters — HBEGF and HBP1 are not haemoglobin."""
    adata = samples["MMRF_1695"].copy()
    qc.add_qc_metrics(adata)
    hb = set(adata.var_names[adata.var["hb"]])
    assert "HBB" in hb and "HBA1" in hb
    assert not {"HBEGF", "HBP1", "HBS1L"} & hb


@requires_data
def test_the_deposit_is_pre_filtered_and_differently_per_cohort(samples):
    """CORRECTS an earlier claim: the depositors' 10,000-UMI cut WAS applied.

    The project plan said their stated QC "demonstrably was not applied to what is
    deposited", reasoning from a cohort-wide UMI average. That average pooled MMRF
    (42% of cells above 10k) with WashU (0% above 10k) and hid a per-cohort truth.
    What the files actually show:

        WU1/WU2   UMI < 10,000   UMI >= 1,000   mt < 20%   genes >= 200
        MMRF      uncensored     UMI >= 1,000   mt < 10%   genes >= 200
        Donor     uncensored     uncensored     mt < 20%   genes >= 200

    This matters well beyond bookkeeping — see the notebook. The censored band is
    enriched 20-40x for GPRC5D, so the WashU cohorts had the antigen-positive tail
    of their own tumours removed before deposit.
    """
    for name, adata in samples.items():
        block = adata.copy()
        qc.add_qc_metrics(block)
        cohort = str(block.obs["cohort"].iat[0])
        counts = block.obs["total_counts"]
        mt = block.obs["pct_counts_mt"]
        genes = block.obs["n_genes_by_counts"]

        assert genes.min() >= 200, f"{name}: gene floor"
        if cohort in {"WU1", "WU2"}:
            assert counts.max() < 10_000, f"{name}: expected the 10k UMI ceiling"
            assert counts.min() >= 1_000, f"{name}: expected the 1k UMI floor"
        if cohort == "MMRF":
            assert counts.max() > 10_000, f"{name}: MMRF is uncensored on UMIs"
            assert mt.max() <= 10.0, f"{name}: MMRF was cut at 10% mt"
        else:
            assert mt.max() <= 20.0, f"{name}: expected the 20% mt ceiling"


@requires_data
def test_our_own_floors_are_a_safety_net_not_the_filter(samples):
    """`min_genes=200` is deliberately below what the deposit already enforces.

    Its job is to catch something that is unambiguously not a cell, not to do the
    filtering — a high floor here would hit the shallow cohorts hardest, which is
    the failure per-cohort thresholds exist to avoid.
    """
    adata = samples["27522_1"].copy()
    qc.add_qc_metrics(adata)
    qc.flag_outliers(adata, group_key=None)
    assert not adata.obs["outlier_min_genes"].any()


@requires_data
def test_sc_best_practices_8pct_mt_cap_would_be_wrong_here(samples):
    """The concrete reason the tutorial's number is not copied.

    A healthy-PBMC 8% cap sits near this marrow's own median, so adopting it would
    discard something close to half the cells for no quality reason.
    """
    adata = samples["BM4"].copy()
    qc.add_qc_metrics(adata)
    assert adata.obs["pct_counts_mt"].median() > 5.0


@requires_data
@requires_r
@pytest.mark.slow
def test_scdblfinder_bridge_runs_and_matches_cell_count(samples):
    adata = samples["BM4"].copy()
    frame = qc.detect_doublets(adata, seed=0)
    assert len(frame) == adata.n_obs
    assert set(frame["class"]) <= {"singlet", "doublet"}
    assert adata.obs["doublet_score"].between(0, 1).all()
    # A plausible rate, not a fixed number — scDblFinder is stochastic.
    rate = (adata.obs["doublet_class"].astype(str) == "doublet").mean()
    assert 0.0 < rate < 0.35


def test_qc_report_columns_do_not_pretend_to_partition():
    """The per-filter counts overlap; only `n_removed` is the union."""
    obs = pd.DataFrame({
        "sample_name": ["A"] * 4,
        "cohort": ["WU1"] * 4,
        "outlier_counts":    [True,  True,  False, False],
        "outlier_genes":     [True,  False, False, False],
        "outlier_top20":     [False, False, False, False],
        "outlier_mt":        [False, False, True,  False],
        "outlier_min_genes": [False, False, False, False],
        "is_doublet":        [False, False, False, False],
        "keep":              [False, False, False, True],
        "n_genes_by_counts": [100.0, 200.0, 300.0, 400.0],
        "total_counts":      [500.0, 600.0, 700.0, 800.0],
        "pct_counts_mt":     [5.0, 6.0, 30.0, 7.0],
    })
    report = qc.qc_report(obs, by="sample_name").iloc[0]
    assert report["n_cells_pre"] == 4
    assert report["n_kept"] == 1
    assert report["n_removed"] == 3
    # Two filters both flag cell 0, so the per-filter counts oversum.
    per_filter = sum(report[c] for c in (
        "n_outlier_counts", "n_outlier_genes", "n_outlier_top20",
        "n_outlier_mt", "n_outlier_min_genes"))
    assert per_filter == 4 > report["n_removed"]
    # Post-QC medians describe the kept cells only.
    assert report["median_genes_post"] == 400.0


@requires_data
@pytest.mark.slow
def test_checkpoints_keep_every_barcode_and_carry_the_cohort_thresholds():
    """Cells are annotated, not deleted — stage 08 must be able to ask what QC cost."""
    import anndata

    from mm_escape import config, io

    directory = config.RESULTS_DIR / "04_qc" / "samples"
    path = directory / "BM4.h5ad"
    if not path.exists():
        pytest.skip("stage 04 has not been run on this machine")

    adata = anndata.read_h5ad(path)
    raw = io.read_sample("BM4")
    assert adata.n_obs == raw.n_obs, "QC must not drop barcodes from the checkpoint"
    for column in ("outlier", "is_doublet", "keep", "doublet_score"):
        assert column in adata.obs
    assert (adata.obs["keep"] == ~(adata.obs["outlier"] | adata.obs["is_doublet"])).all()
    # Written back by the cohort pass, so the file carries the shipped thresholds.
    assert adata.uns["qc_threshold_scope"] == "cohort"
    assert set(adata.uns["qc_thresholds"]["group"]) >= {"Donor"}
