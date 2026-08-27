# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 06
# Step         : s06a_cluster_de_top_genes.sh
# What it does : cluster-level DE marker table
# Writes       : results/06_annotation/cluster_de_top_genes.csv
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T05:01:37Z
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

# --- add gc + RSS helper to the imports cell ---
s=s.replace(
"""from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path""",
"""from __future__ import annotations

import gc
import json
import resource
import sys
import warnings
from pathlib import Path""",1)

s=s.replace(
'''sc.settings.verbosity = 1
CLUSTER_KEY = "leiden"''',
'''sc.settings.verbosity = 1
CLUSTER_KEY = "leiden"


def rss(tag: str) -> None:
    """Peak RSS so far. This stage loads the same object stage 05 peaked ~20 GB on,
    and the box has 30 GB — memory is a real constraint here, not a footnote."""
    gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
    print(f"[mem] peak RSS {gb:5.1f} GB  after {tag}", flush=True)''',1)

s=s.replace('print(adata.shape)\nprint(adata.obs[CLUSTER_KEY].value_counts().sort_index().to_string())',
            'print(adata.shape)\nprint(adata.obs[CLUSTER_KEY].value_counts().sort_index().to_string())\nrss("load")',1)

# --- move the SingleR pseudobulk + CellTypist input up, then free the counts layer ---
old_ct = '''import celltypist                                # noqa: E402
from celltypist import models                    # noqa: E402

ct_input = adata.copy()
ct_input.X = ct_input.layers["counts"].copy()
sc.pp.normalize_total(ct_input, target_sum=1e4)
sc.pp.log1p(ct_input)

models.download_models(model=["Immune_All_Low.pkl", "Immune_All_High.pkl"], force_update=False)'''
new_ct = '''import anndata as ad                             # noqa: E402
import celltypist                                # noqa: E402
from celltypist import models                    # noqa: E402

# Everything that needs `layers["counts"]` happens here, then the layer is dropped.
# A full `adata.copy()` here is what killed the first run: it duplicates a
# 172,940 x 32,991 object on a 30 GB box that stage 05 already peaked ~20 GB on.
# Build a lean object that carries only what CellTypist reads instead.
pseudo = sc.get.aggregate(adata, by=CLUSTER_KEY, func="mean", layer="counts")
pseudo_df = pd.DataFrame(
    pseudo.layers["mean"], index=pseudo.obs_names, columns=pseudo.var_names
).T
print("cluster pseudobulk for SingleR:", pseudo_df.shape)

ct_input = ad.AnnData(
    X=adata.layers["counts"].copy(),
    obs=pd.DataFrame(index=adata.obs_names),
    var=pd.DataFrame(index=adata.var_names),
)
del adata.layers["counts"]
gc.collect()
rss("counts layer freed")

sc.pp.normalize_total(ct_input, target_sum=1e4)
sc.pp.log1p(ct_input)

models.download_models(model=["Immune_All_Low.pkl", "Immune_All_High.pkl"], force_update=False)'''
assert old_ct in s
s=s.replace(old_ct,new_ct,1)

s=s.replace('''ct_lab = ct_res.predicted_labels
adata.obs["celltypist_raw"] = ct_lab["majority_voting"].to_numpy()
adata.obs["celltypist_conf"] = ct_res.probability_matrix.max(axis=1).to_numpy()
del ct_input''',
'''ct_lab = ct_res.predicted_labels
adata.obs["celltypist_raw"] = ct_lab["majority_voting"].to_numpy()
adata.obs["celltypist_conf"] = ct_res.probability_matrix.max(axis=1).to_numpy()
del ct_input, ct_res
gc.collect()
rss("celltypist")''',1)

# --- SingleR cell no longer recomputes the pseudobulk ---
old_sr = '''pseudo = sc.get.aggregate(adata, by=CLUSTER_KEY, func="mean", layer="counts")
pseudo_df = pd.DataFrame(pseudo.layers["mean"], index=pseudo.obs_names, columns=pseudo.var_names).T
print("cluster pseudobulk for SingleR:", pseudo_df.shape)

ro.r('''
assert old_sr in s
s=s.replace(old_sr,"""ro.r('""",1)

# --- marker scoring is where the DE pass follows; report memory around both ---
s=s.replace('print("marker scores:", score_cols)','print("marker scores:", score_cols)\nrss("marker scoring")',1)
s=s.replace('de_top.to_csv(OUT / "cluster_de_top_genes.csv")',
            'de_top.to_csv(OUT / "cluster_de_top_genes.csv")\nrss("wilcoxon DE")',1)
s=s.replace('print(f"\\nveto line: {config.MARKER_COVERAGE_MIN}")',
            'print(f"\\nveto line: {config.MARKER_COVERAGE_MIN}")\nrss("marker coverage")',1)
open(p,'w').write(s)
print("notebook restructured for memory")
PY
grep -n "rss(\|del adata.layers\|ad.AnnData(\|adata.copy()" notebooks/06_annotation.py
