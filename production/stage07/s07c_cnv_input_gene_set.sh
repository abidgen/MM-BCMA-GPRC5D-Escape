# ============================================================================
# RECOVERED PRODUCTION DRIVER  --  verbatim historical source, NOT re-executed
#
# Stage        : 07
# Step         : s07c_cnv_input_gene_set.sh
# What it does : CNV input gene set and exclusions (antigens excluded by construction)
# Writes       : results/07_malignant_plasma/cnv_negative_result/cnv_input_gene_set.tsv.gz; cnv_excluded_genes.tsv.gz
#
# PROVENANCE
#   Recovered from : Claude Code session transcript
#                    ~/.claude/projects/-media-wrath-CART-mm-dual-antigen/77d8225f-3278-4204-8b06-7f58a1d38646.jsonl
#   Original form  : inline shell/python block
#   Original path  : (inline heredoc executed directly from the shell)
#   Executed (UTC) : 2026-08-25T11:52:51Z
#   Recovered on   : 2026-08-26 (pre-Stage-12 provenance repair)
#
# STATUS: PRODUCER_RECOVERED_VERBATIM.  This is the code that actually produced
# the frozen artifacts listed above.  It has NOT been re-run during recovery and
# no frozen artifact was regenerated.  Re-running it OVERWRITES frozen results --
# see production/README.md before executing anything in this tree.
# ============================================================================
cd /media/wrath/CART_mm_dual_antigen
/home/abid/miniforge3/envs/mm-core/bin/python << 'EOF' 2>&1 | grep -v Warning
import pandas as pd, numpy as np, re
pos = pd.read_csv('raw/gtf/gene_space_positions.tsv.gz', sep='\t')
print(f"gene space with positions: {len(pos)}")

# --- IG loci: symbol-based (V/D/J/C) + coordinate verification -------------------
IG_RX = re.compile(r'^IG[HKL](V|D|J|C)')
CONST = {'IGKC','IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','IGHM','IGHD','IGHE',
         'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7'}
sym = pos['symbol'].astype(str)
is_ig = sym.str.match(IG_RX) | sym.isin(CONST)

# canonical GRCh38 locus windows
LOCI = {'IGH (14q32)': ('14', 105_580_000, 106_880_000),
        'IGK (2p11)':  ('2',  88_780_000,  90_280_000),
        'IGL (22q11)': ('22', 22_020_000,  22_930_000)}
in_locus = pd.Series(False, index=pos.index)
for name,(c,s,e) in LOCI.items():
    m = (pos.chromosome.astype(str)==c) & (pos.start>=s) & (pos.end<=e)
    in_locus |= m
    print(f"  {name}: {int(m.sum())} genes in window | IG-symbol among them {int((m & is_ig).sum())}")

excl = is_ig | in_locus
print(f"\nIG genes by symbol pattern : {int(is_ig.sum())}")
print(f"IG genes inside loci windows: {int(in_locus.sum())}")
print(f"UNION excluded              : {int(excl.sum())}")
print("\nper-locus breakdown of the union (by chromosome):")
print(pos[excl].groupby(pos.loc[excl,'chromosome'].astype(str)).size().sort_values(ascending=False).head(8).to_string())
print("\nIG-symbol genes NOT in a canonical window (orphons/other chromosomes):")
orph = pos[is_ig & ~in_locus]
print(f"  {len(orph)} genes on chromosomes {sorted(orph.chromosome.astype(str).unique())}")

# --- antigens ---------------------------------------------------------------
ANT = ['TNFRSF17','GPRC5D']
is_ant = sym.isin(ANT)
print(f"\nantigen genes excluded: {sorted(sym[is_ant])}")

CANON = [str(x) for x in range(1,23)] + ['X']
keep = (~excl) & (~is_ant) & pos.chromosome.astype(str).isin(CANON)
cnv = pos[keep].copy()
print(f"\n=== CNV INPUT GENE SET ===")
print(f"  start 32,991 -> IG {int(excl.sum())} -> antigen {int(is_ant.sum())} "
      f"-> non-canonical/MT/Y {int(((~excl)&(~is_ant)&~pos.chromosome.astype(str).isin(CANON)).sum())}")
print(f"  FINAL: {len(cnv)} genes")
cov = cnv.groupby(cnv.chromosome.astype(str)).size().reindex(CANON)
print(f"\n  per-chromosome coverage after removal (min {int(cov.min())} on chr{cov.idxmin()}):")
print("   " + "  ".join(f"{c}:{int(n)}" for c,n in cov.items()))
print(f"  chromosomes with <100 genes: {list(cov[cov<100].index) or 'none'}")
cnv.sort_values(['chromosome','start']).to_csv('results/07_malignant_plasma/cnv_input_gene_set.tsv.gz',
                                               sep='\t', index=False, compression='gzip')
pos[excl | is_ant].to_csv('results/07_malignant_plasma/cnv_excluded_genes.tsv.gz',
                          sep='\t', index=False, compression='gzip')
print("\n  wrote cnv_input_gene_set.tsv.gz + cnv_excluded_genes.tsv.gz")
EOF
