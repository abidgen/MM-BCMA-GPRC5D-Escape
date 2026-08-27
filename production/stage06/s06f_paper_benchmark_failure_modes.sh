# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 06
# Step         : s06f_paper_benchmark_failure_modes.sh
# What it does : paper-annotation benchmark, failure modes
# Writes       : results/06_annotation/paper_annotation_benchmark/*
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T10:49:53Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
/home/abid/miniforge3/envs/mm-annotation/bin/python << 'EOF' 2>&1 | grep -v Warning | tail -50
import scanpy as sc, numpy as np, pandas as pd
from pathlib import Path
OUT=Path('results/06_annotation/paper_annotation_benchmark')
pd.set_option('display.width',250)
a=sc.read_h5ad('results/06_annotation_c2d_accepted/annotated.h5ad')
cl=a.obs['leiden'].astype(str)
gs=['CD3E','CD3D','CD8A','CD8B','CD7','CD4','IL7R','TRAC','TRDC','NKG7','GNLY',
    'KLRD1','KLRF1','NCAM1','FCGR3A','HBB','HBA1','AHSP','GYPA','ALAS2','CD34','SPINK2','SOX4',
    'CD14','LYZ','FCER1A','CLEC10A','CLEC4C','LILRA4','MZB1','SDC1','IGHG1','XBP1','MS4A1','CD79A']
gs=[g for g in gs if g in a.var_names]
X=a[:,gs].X; X=np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
det=pd.DataFrame(X>0,columns=gs).groupby(cl.values).mean()*100

comp=pd.read_csv(OUT/'cluster_paper_comparison.csv',index_col=0)
KEY=[0,3,4,7,10,12,21,22,23,24,25,26]
k=comp.loc[KEY].copy()
k.to_csv(OUT/'key_cluster_comparison.csv')
print("=== KEY CLUSTERS ==="); print(k[['n_cells','paper_fine','paper_broad','c2d_manual','c2d_final','agrees_broad']].to_string())
print()
print("=== FM1: NKG7/GNLY shared by NK and cytotoxic T (clusters 3,12,23,19) ===")
print(det.loc[['3','12','23','19'],['CD3E','CD3D','CD8A','CD8B','CD7','TRAC','TRDC','NKG7','GNLY','KLRF1','NCAM1']].round(1).to_string())
print()
print("=== FM2: globin/AHSP ambient vs true erythroid (8,24 disputed; 27 true) ===")
print(det.loc[['8','24','27','1'],['HBB','HBA1','AHSP','GYPA','ALAS2','CD34','SPINK2','SOX4']].round(1).to_string())
print()
print("=== FM3/FM4/FM5: CD14/LYZ insufficiency; DC vs B; broad innate (22,25,26,10,21) ===")
print(det.loc[['22','25','26','10','21'],['CD14','LYZ','FCGR3A','FCER1A','CLEC10A','CLEC4C','LILRA4','MS4A1','CD79A']].round(1).to_string())
print()
print("=== FM6: MZB1/XBP1 overlap plasma vs progenitor/B (0,4,7 true plasma; 24,21,26) ===")
print(det.loc[['0','4','7','24','21','26'],['MZB1','XBP1','SDC1','IGHG1','CD34','SPINK2','MS4A1']].round(1).to_string())
det.round(2).to_csv(OUT/'failure_mode_marker_detection.csv')
EOF
