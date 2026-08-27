# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 06
# Step         : s06b_lineage_contradiction.sh
# What it does : per-method lineage-contradiction tables
# Writes       : results/06_annotation/lineage_contradiction_*.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T05:27:12Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
python3 - <<'PY'
p='notebooks/06_annotation.py'
s=open(p).read()

# --- load cell: reuse v1's automated predictions when inputs are unchanged ---
old='''adata = sc.read_h5ad(IN_H5AD)
print(adata.shape)
print(adata.obs[CLUSTER_KEY].value_counts().sort_index().to_string())
rss("load")'''
new='''# Stage 06 ran twice. v1 is preserved at `results/06_annotation_v1/` and was NOT
# accepted — see CLAUDE.md and docs/decisions-archive.md. v2 revises the manual
# reference and adds the lineage-exclusivity veto; the acceptance thresholds are
# unchanged.
#
# CellTypist and SingleR are reused from v1 rather than recomputed. That is legitimate
# here and only here: their inputs (the stage-05 expression matrix), their models and
# references, and the prediction procedure are all untouched by this revision, which
# changes only the MANUAL panel and the validation framework. The barcode identity is
# asserted below, so a silently different input fails loudly instead of being reused.
V1 = REPO / "results" / "06_annotation_v1" / "annotated.h5ad"
REUSED = ("cell_type_celltypist", "cell_type_singler_nov", "cell_type_singler_hpca",
          "celltypist_raw", "celltypist_percell", "celltypist_conf",
          "singler_nov_raw", "singler_hpca_raw")

if V1.exists():
    adata = sc.read_h5ad(V1)
    missing = [c for c in REUSED if c not in adata.obs]
    assert not missing, f"v1 predictions incomplete: {missing}"
    ref_obs = sc.read_h5ad(IN_H5AD, backed="r").obs_names
    assert list(adata.obs_names) == list(ref_obs), "v1 barcodes differ from stage 05 — do not reuse"
    # Drop everything v2 recomputes, so no v1 result can leak into a v2 table.
    for c in ["cell_type_manual", "cell_type", "cell_type_fine", "annotation_source",
              "annotation_conf", *[c for c in adata.obs.columns
                                   if c.startswith(("score_", "program_"))]]:
        if c in adata.obs:
            del adata.obs[c]
    REUSE = True
    print(f"reusing v1 automated predictions for {len(REUSED)} columns")
else:
    adata = sc.read_h5ad(IN_H5AD)
    REUSE = False

print(adata.shape)
print(adata.obs[CLUSTER_KEY].value_counts().sort_index().to_string())
rss("load")'''
assert old in s
s=s.replace(old,new,1)

# --- guard the CellTypist block ---
old_ct='''import anndata as ad                             # noqa: E402
import celltypist                                # noqa: E402
from celltypist import models                    # noqa: E402
'''
new_ct='''import anndata as ad                             # noqa: E402
'''
assert old_ct in s
s=s.replace(old_ct,new_ct,1)

start = s.index('# Everything that needs `layers["counts"]` happens here')
end = s.index('rss("celltypist")') + len('rss("celltypist")')
block = s[start:end]
s = s[:start] + '''if not REUSE:
    import celltypist                                        # noqa: E402
    from celltypist import models                            # noqa: E402
''' + "\n".join("    " + ln if ln.strip() else ln for ln in block.split("\n")) + '''
else:
    pseudo_df = None      # SingleR is reused too; no pseudobulk needed
    if "counts" in adata.layers:
        del adata.layers["counts"]
        gc.collect()
    rss("reused celltypist")
''' + s[end:]

# --- guard the SingleR block ---
sr_start = s.index('# rpy2 resolves R via $R_HOME')
sr_end = s.index('# %% [markdown]\n# ## Marker coverage')
sr_block = s[sr_start:sr_end]
s = s[:sr_start] + '''if not REUSE:
''' + "\n".join("    " + ln if ln.strip() else ln for ln in sr_block.rstrip().split("\n")) + '''
else:
    singler_out = {"NovershternHematopoieticData": "reused from v1",
                   "HumanPrimaryCellAtlasData": "reused from v1"}
    for short in ("singler_nov", "singler_hpca"):
        print(f"{short}: reused")
        print(adata.obs[f"cell_type_{short}"].value_counts().to_string())

''' + s[sr_end:]

# --- add the exclusivity computation before the concordance section ---
old_cov='''print(f"\\nveto line: {config.MARKER_COVERAGE_MIN}")
rss("marker coverage")'''
new_cov='''print(f"\\nveto line: {config.MARKER_COVERAGE_MIN}")
rss("marker coverage")

# %% [markdown]
# ### Lineage exclusivity — the second veto (v2)
#
# Coverage asks whether cells labelled X express X's markers, which is
# precision-like and therefore blind to a class that has *swallowed* another lineage.
# This asks the complementary question: do those same cells carry strong positive
# evidence for a lineage incompatible with X? Detection-based, so dropout can only
# hide a contradiction and never invent one.

# %%
contradiction = {m: ann.contradiction_rate(adata, k) for m, k in METHODS.items()}
contra_wide = pd.DataFrame({m: t["contradiction_rate"] for m, t in contradiction.items()})
contra_wide.to_csv(OUT / "lineage_contradiction_by_method.csv")
for m, t in contradiction.items():
    t.to_csv(OUT / f"lineage_contradiction_{m}.csv")
print(contra_wide.round(3).to_string())
print(f"\\nveto line: {config.CONTRADICTION_MAX_RATE}")
rss("lineage exclusivity")'''
assert old_cov in s
s=s.replace(old_cov,new_cov,1)

s=s.replace('decision = ann.decide_per_class(conc, coverage)',
            'decision = ann.decide_per_class(conc, coverage, contradiction)',1)
s=s.replace('''print(decision[["chosen_method", "reason", "f1_threshold", "f1_chosen",
                "coverage_chosen", "vetoed_by_coverage"]].to_string())''',
            '''print(decision[["chosen_method", "reason", "f1_threshold", "f1_chosen",
                "coverage_chosen", "contradiction_chosen",
                "vetoed_by_coverage", "vetoed_by_contradiction"]].to_string())''',1)
s=s.replace('''    "decision": decision["chosen_method"].to_dict(),''',
            '''    "decision": decision["chosen_method"].to_dict(),
    "contradiction_max": config.CONTRADICTION_MAX_RATE,
    "contradiction_min_genes": config.CONTRADICTION_MIN_GENES,
    "reused_v1_predictions": bool(REUSE),''',1)
s=s.replace('OUT = REPO / "results" / "06_annotation"','OUT = REPO / "results" / "06_annotation"',1)
open(p,'w').write(s)
print("v2 notebook written")
PY
/home/abid/miniforge3/envs/mm-annotation/bin/python -c "import ast;ast.parse(open('notebooks/06_annotation.py').read());print('syntax OK')"
