# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 06
# Step         : s06e_paper_benchmark_markers.sh
# What it does : paper-annotation benchmark, marker detection
# Writes       : results/06_annotation/paper_annotation_benchmark/*
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T10:49:25Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
mkdir -p results/06_annotation/paper_annotation_benchmark
/home/abid/miniforge3/envs/mm-annotation/bin/python << 'EOF' 2>&1 | grep -v Warning | tail -45
import scanpy as sc, numpy as np, pandas as pd
from pathlib import Path
OUT = Path('results/06_annotation/paper_annotation_benchmark')

# Paper markers, verbatim from Materials and Methods, "scRNA-seq cell type annotation"
# (Cancer Research 2023;83:1214-1233, PMC10102848). Cell types assigned by MANUAL
# review of marker expression -- no SingleR, no reference atlas.
PAPER = {
 "B cells":       ["CD79A","CD79B","MS4A1"],
 "CD8+ T cells":  ["CD8A","CD8B","CD7","CD3E"],
 "CD4+ T cells":  ["CD4","IL7R","CD7","CD3E"],
 "NK cells":      ["NKG7","GNLY"],
 "Plasma cells":  ["MZB1","SDC1","IGHG1"],
 "Macrophages":   ["FCGR3A"],
 "Monocytes":     ["CD14","LYZ"],
 "Dendritic cells":["FCER1A","CLEC10A"],
 "Erythrocytes":  ["AHSP","HBA1","HBA2","HBB"],   # paper writes "AHSP1, HBA, HBB"
 "pDC":           ["CLEC4C","LILRA4"],
 "Neutrophils":   ["AZU1","MPO","ELANE"],
}
COLLAPSE = {"B cells":"Bcell","CD8+ T cells":"Tcell","CD4+ T cells":"Tcell",
 "NK cells":"NK","Plasma cells":"PlasmaCell","Macrophages":"Myeloid",
 "Monocytes":"Myeloid","Dendritic cells":"Myeloid","pDC":"Myeloid",
 "Erythrocytes":"Erythroid","Neutrophils":"Myeloid"}

a = sc.read_h5ad('results/06_annotation_c2d_accepted/annotated.h5ad')
cl = a.obs['leiden'].astype(str)
genes = sorted({g for v in PAPER.values() for g in v})
present = [g for g in genes if g in a.var_names]
X = a[:, present].X
X = np.asarray(X.todense()) if hasattr(X,'todense') else np.asarray(X)
det = pd.DataFrame(X>0, columns=present).groupby(cl.values).mean()

rows=[]
for ct, gs in PAPER.items():
    hit=[g for g in gs if g in present]
    rows.append({"paper_cell_type":ct,"reported_markers":"; ".join(gs),
      "markers_in_gene_space":"; ".join(hit),
      "missing":"; ".join(sorted(set(gs)-set(hit))) or "none",
      "collapses_to":COLLAPSE[ct],
      "source":"Methods, 'scRNA-seq cell type annotation' (PMC10102848)"})
pd.DataFrame(rows).to_csv(OUT/'paper_marker_mapping.csv', index=False)

# Paper-style call: argmax of mean marker detection -- their "manual review" rule.
score = pd.DataFrame({ct: det[[g for g in gs if g in present]].mean(axis=1)
                      for ct, gs in PAPER.items()})
paper_fine = score.idxmax(axis=1)
paper_broad = paper_fine.map(COLLAPSE)

o = a.obs
def mode(col, c):
    s = o.loc[(cl==c).values, col].astype(str)
    return s.mode().iat[0] if len(s) else '-'
comp=[]
for c in det.index:
    comp.append({"leiden":int(c), "n_cells":int((cl==c).sum()),
      "paper_fine":paper_fine[c], "paper_broad":paper_broad[c],
      "paper_top_score":round(float(score.loc[c].max()),3),
      "paper_runner_up":score.loc[c].nlargest(2).index[-1],
      "paper_runner_up_score":round(float(score.loc[c].nlargest(2).iloc[-1]),3),
      "c2d_manual":mode('cell_type_manual',c), "celltypist":mode('cell_type_celltypist',c),
      "singler_nov":mode('cell_type_singler_nov',c), "singler_hpca":mode('cell_type_singler_hpca',c),
      "c2d_final":mode('cell_type',c)})
comp=pd.DataFrame(comp).sort_values('leiden')
comp["agrees_broad"]=comp.paper_broad==comp.c2d_final
comp.to_csv(OUT/'cluster_paper_comparison.csv', index=False)
det.round(4).to_csv(OUT/'paper_marker_detection_by_cluster.csv')

print("=== paper markers missing from gene space ===")
print(sorted(set(genes)-set(present)) or "none")
print()
print("=== all 30 clusters ===")
print(comp[['leiden','n_cells','paper_fine','paper_broad','c2d_final','agrees_broad']].to_string(index=False))
print()
print("agreement:", int(comp.agrees_broad.sum()), "/", len(comp), "clusters;",
      f"{100*comp.loc[comp.agrees_broad,'n_cells'].sum()/comp.n_cells.sum():.1f}% of cells")
EOF
