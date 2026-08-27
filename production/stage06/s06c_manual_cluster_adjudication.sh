# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 06
# Step         : s06c_manual_cluster_adjudication.sh
# What it does : manual cluster adjudication
# Writes       : results/06_annotation/manual_cluster_adjudication.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T05:46:58Z
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
s=s.replace('V1 = REPO / "results" / "06_annotation_v1" / "annotated.h5ad"',
            'V_PREV = REPO / "results" / "06_annotation_v2" / "annotated.h5ad"',1)
s=s.replace('if V1.exists():\n    adata = sc.read_h5ad(V1)','if V_PREV.exists():\n    adata = sc.read_h5ad(V_PREV)',1)
s=s.replace('assert not missing, f"v1 predictions incomplete: {missing}"',
            'assert not missing, f"previous predictions incomplete: {missing}"',1)
s=s.replace('assert list(adata.obs_names) == list(ref_obs), "v1 barcodes differ from stage 05 — do not reuse"',
            'assert list(adata.obs_names) == list(ref_obs), "barcodes differ from stage 05 — do not reuse"',1)
s=s.replace('print(f"reusing v1 automated predictions for {len(REUSED)} columns")',
            'print(f"reusing v2 automated predictions for {len(REUSED)} columns")',1)

old='''manual_tbl = ann.manual_labels_from_clusters(adata, CLUSTER_KEY, out_key="cell_type_manual")
manual_tbl.to_csv(OUT / "manual_cluster_assignment.csv")
print(manual_tbl[["winner", "best", "second", "margin"]].to_string())'''
new='''# v3: identity is adjudicated on DETECTION FRACTIONS, not on cross-panel score_genes
# argmax. Those module scores subtract a control set drawn from each gene's own
# expression bin, so each panel carries its own baseline offset — measured here as
# 0.2036 favouring NK over T, larger than the 0.068/0.196 margins by which two T-cell
# clusters were called NK in v1/v2. The scores above stay as descriptive within-program
# quantities and no longer decide anything.
detection = ann.marker_detection_by_cluster(adata, CLUSTER_KEY)
detection.to_csv(OUT / "marker_detection_by_cluster.csv")

manual_tbl = ann.adjudicate_clusters(adata, CLUSTER_KEY, out_key="cell_type_manual")
manual_tbl.to_csv(OUT / "manual_cluster_adjudication.csv")
print(manual_tbl[["winner", "reason", "lead_positive", "runner_up_positive",
                  "supported", "survivors"]].to_string())'''
assert old in s
s=s.replace(old,new,1)
s=s.replace('''    "reused_v1_predictions": bool(REUSE),''',
            '''    "reused_v2_predictions": bool(REUSE),
    "manual_adjudication": {"detect_min": config.MANUAL_MARKER_DETECT_MIN,
                            "positive_min": config.MANUAL_POSITIVE_MIN,
                            "margin": config.MANUAL_DECISION_MARGIN},''',1)
open(p,'w').write(s)
print("notebook wired for v3")
PY
/home/abid/miniforge3/envs/mm-annotation/bin/python -c "import ast;ast.parse(open('notebooks/06_annotation.py').read());print('syntax OK')"
/home/abid/miniforge3/envs/mm-annotation/bin/python -u notebooks/06_annotation.py > /tmp/claude-1000/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646/scratchpad/s06/run06_v3.log 2>&1
