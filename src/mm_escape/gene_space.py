"""
Cross-reference gene-space harmonization for GSE223060.

THE PROBLEM (confirmed ground truth — do not re-derive)
-------------------------------------------------------
The 62 samples were processed against three different Cell Ranger references,
distinguishable by row count in `genes.tsv`:

    33538 genes -> 37 samples   refdata-cellranger-GRCh38-3.0.0 (Ensembl 93)
    33694 genes -> 24 samples   refdata-cellranger-GRCh38-1.2.0 (Ensembl 84)
    22185 rows   ->  1 sample    56203_1 only — a TRUNCATED 33694 file, repaired
                                on read by io.read_sample (see config.TRUNCATED_GENE_FILES)

The two retained references use different HGNC symbol vintages, so intersecting on
gene symbols silently drops genes present in BOTH builds — including NSD2 (WHSC1),
without which t(4;14) becomes uncallable and stage 10 loses its highest-risk TC class.

Worse, the deposited files went through Seurat, which applied `gsub("_", "-")` and
R's `make.unique`. The resulting `.N` suffixes encode ROW ORDER, not biology, so a
symbol join can pair the WRONG gene rather than merely miss one:

    33538 build:  TBCE -> ENSG00000285053
    33694 build:  TBCE -> ENSG00000116957     <- a different annotation entry

THE FIX
-------
Join on Ensembl gene ID. The deposit has no ID column — `genes.tsv` is a single column
of symbols and contains zero `ENSG` strings across all 62 samples — but there are only
three distinct gene files in the cohort (checksum-verified, byte-identical within
group), each a positional dump of a public reference. So the IDs are reconstructible:

    1. Ensembl GTF (release 93 / 84), `feature == "gene"` rows in GTF order,
       first occurrence per gene_id, filtered to 10x's mkgtf biotype whitelist
       -> exactly 33538 / 33694 rows.
    2. Apply Seurat's transforms in order: gsub("_", "-"), then R make.unique.
    3. ASSERT the result equals the deposited column position-for-position.

Step 3 is what makes this self-certifying rather than a guess: a wrong release or
biotype filter changes the row count or the order, and the assertion fails loudly.
Verified 2026-08-21: 0 mismatches / 33538 rows and 0 mismatches / 33694 rows.

PAYOFF
------
    raw symbols                 22,164 genes
    symbols + 4-gene alias map  22,168 genes
    Ensembl IDs                 32,991 genes   (+10,827, ~49% more gene space)

11,140 intersected IDs carry a different symbol in each build and were invisible to a
symbol join.

INDEX CONVENTION
----------------
`ensembl_id` is `var_names` THROUGH THE MERGE ONLY — that is where identity is
load-bearing. Once the object is a single harmonized matrix the mis-pairing risk is
gone, so `to_canonical_symbols()` switches `var_names` to the Ensembl-93 symbol and
retains the IDs in `.var`. Everything downstream (score_genes panels, dotplots,
celltypist, decoupler, liana) is symbol-native.

NEVER call `var_names_make_unique()` on these objects. It would re-apply exactly the
positional mangling this module exists to undo. The 9 symbols that still collide after
the ID intersection are suffixed `SYMBOL__ENSGxxxxxxxxxxx` instead.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from . import config

__all__ = [
    "load_gene_map",
    "detect_build",
    "attach_ensembl_ids",
    "intersect_gene_space",
    "to_canonical_symbols",
    "assert_required_genes",
    "rebuild_gene_map_from_gtf",
    "GeneSpaceError",
]


class GeneSpaceError(RuntimeError):
    """Raised when the gene space is not what the pipeline requires.

    Always carries the specific offending gene/sample names — a silent partial
    marker set is the failure mode this whole module exists to prevent.
    """


_ENSEMBL_ID_RE = re.compile(r"^ENSG\d{11}$")
_VERSION_SUFFIX_RE = re.compile(r"\.\d+$")

_GTF_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_GTF_GENE_NAME = re.compile(r'gene_name "([^"]+)"')
_GTF_GENE_BIOTYPE = re.compile(r'gene_biotype "([^"]+)"')


# ---------------------------------------------------------------------------
# The two Seurat transforms the depositor applied
# ---------------------------------------------------------------------------

def _seurat_sanitize(symbol: str) -> str:
    """Seurat replaces underscores with hyphens in feature names.

    This is why the deposit spells Ensembl's `RP11-442N24__B` as `RP11-442N24--B`.
    Applies to 15 rows in the 33694 build, 0 in the 33538 build.
    """
    return symbol.replace("_", "-")


def _make_unique(names: Sequence[str]) -> list[str]:
    """Reimplementation of R's `make.unique` (dot-suffix disambiguation).

    First occurrence keeps the bare name; subsequent ones get `.1`, `.2`, ...,
    skipping any candidate already taken. This is NOT the same as
    `anndata.AnnData.var_names_make_unique` and the two must not be mixed —
    reproducing the depositor's exact output is the point.
    """
    counter: dict[str, int] = defaultdict(int)
    taken: set[str] = set()
    out: list[str] = []
    for name in names:
        if name not in taken:
            out.append(name)
            taken.add(name)
            counter[name] = 0
            continue
        k = counter[name]
        while True:
            k += 1
            candidate = f"{name}.{k}"
            if candidate not in taken:
                break
        counter[name] = k
        taken.add(candidate)
        out.append(candidate)
    return out


# ---------------------------------------------------------------------------
# Loading the committed map
# ---------------------------------------------------------------------------

def load_gene_map(build: int, gene_space_dir: Path | None = None) -> pd.DataFrame:
    """Load the committed `deposited_symbol -> ensembl_id` map for one build.

    Returns a DataFrame indexed by `row_index` with columns
    `deposited_symbol`, `ensembl_id`, ordered exactly as the deposited `genes.tsv`.

    The map is committed at `resources/gene_space/` (~1 MB gzipped) precisely so the
    41-44 MB GTFs never need re-downloading. Regenerate only via
    `rebuild_gene_map_from_gtf`, which re-verifies against the real files.
    """
    if build not in config.BUILDS:
        raise GeneSpaceError(
            f"Unknown reference build with {build} genes. Known: "
            f"{sorted(config.BUILDS)}. A new row count means a new reference — "
            f"reconstruct and verify it before using the sample."
        )

    directory = gene_space_dir or config.GENE_SPACE_DIR
    path = directory / f"genes_{build}_ensembl.tsv.gz"
    if not path.exists():
        raise GeneSpaceError(
            f"Missing gene map {path}. Rebuild with rebuild_gene_map_from_gtf("
            f"gtf_path, build={build}) using {config.BUILDS[build]['gtf_url']}"
        )

    frame = pd.read_csv(path, sep="\t", index_col="row_index")
    if len(frame) != build:
        raise GeneSpaceError(
            f"Gene map {path} has {len(frame)} rows, expected {build}. "
            f"The map is corrupt — regenerate it."
        )
    return frame


def detect_build(var_names: Sequence[str] | pd.Index) -> int:
    """Identify which reference an object came from, by gene count.

    Row count is the discriminator confirmed against all 62 samples; there are only
    three distinct gene files in the cohort and they differ in length.
    """
    n = len(var_names)
    if n in config.BUILDS:
        return n
    truncated = {v["deposited_rows"] for v in config.TRUNCATED_GENE_FILES.values()}
    hint = (
        " That is the row count of a known TRUNCATED deposit, which io.read_sample "
        "repairs on read — load the sample through it rather than by hand."
        if n in truncated else " Was this object already subset or intersected?"
    )
    raise GeneSpaceError(
        f"{n} genes matches no known reference build ({sorted(config.BUILDS)})." + hint
    )


# ---------------------------------------------------------------------------
# Attaching IDs to an AnnData
# ---------------------------------------------------------------------------

def attach_ensembl_ids(adata, gene_space_dir: Path | None = None, copy: bool = False):
    """Set `var_names` to Ensembl IDs, verifying the deposited symbols first.

    Call this on each per-sample AnnData BEFORE any concatenation. The verification
    is position-for-position against the committed map: if the object's `var_names`
    are not the exact deposited column in the exact deposited order, this raises
    rather than producing a plausible-looking wrong join.

    Adds to `.var`:
        deposited_symbol  — the symbol as it appeared in genes.tsv (Seurat-mangled)
        ensembl_id        — the recovered stable identifier

    Note `copy=False` mutates in place, which is what you want when looping over 61
    samples; pass `copy=True` in a notebook when inspecting.
    """
    if copy:
        adata = adata.copy()

    build = detect_build(adata.var_names)
    gene_map = load_gene_map(build, gene_space_dir=gene_space_dir)

    deposited = list(adata.var_names)
    expected = list(gene_map["deposited_symbol"])
    if deposited != expected:
        mismatches = [
            (i, got, want)
            for i, (got, want) in enumerate(zip(deposited, expected))
            if got != want
        ]
        head = mismatches[:5]
        raise GeneSpaceError(
            f"var_names do not match the {build}-gene reference position-for-position "
            f"({len(mismatches)} mismatches). First: "
            + "; ".join(f"row {i}: got {got!r}, expected {want!r}" for i, got, want in head)
            + ". The object was reordered, renamed, or var_names_make_unique() was "
            "called on it — any of which breaks the positional join."
        )

    adata.var["deposited_symbol"] = deposited
    adata.var["ensembl_id"] = list(gene_map["ensembl_id"])
    adata.var_names = pd.Index(gene_map["ensembl_id"].to_numpy(), name="ensembl_id")
    return adata


# ---------------------------------------------------------------------------
# The intersection
# ---------------------------------------------------------------------------

def intersect_gene_space(adatas: Iterable, verbose: bool = True) -> list:
    """Subset every object to the Ensembl IDs shared by all of them.

    INTERSECT, NEVER UNION. A union merge would make ~11k genes structurally zero in
    whole sample cohorts — indistinguishable downstream from a true biological zero,
    which is exactly the quantity this project measures.

    Objects must already have Ensembl IDs as `var_names` (see `attach_ensembl_ids`).
    Order of the returned list matches the input; each object is subset to the same
    ID order, so `anndata.concat(..., join="inner")` is then a no-op on the var axis.
    """
    objects = list(adatas)
    if not objects:
        raise GeneSpaceError("No objects to intersect.")

    for i, adata in enumerate(objects):
        bad = [v for v in adata.var_names[:50] if not _ENSEMBL_ID_RE.match(str(v))]
        if bad:
            raise GeneSpaceError(
                f"Object {i} is not keyed on Ensembl IDs (saw {bad[:3]}). "
                f"Call attach_ensembl_ids() on every object before intersecting — "
                f"intersecting on symbols is the bug this module exists to prevent."
            )

    shared: set[str] = set(objects[0].var_names)
    for adata in objects[1:]:
        shared &= set(adata.var_names)
    if not shared:
        raise GeneSpaceError("Gene-space intersection is empty.")

    # Deterministic order, taken from the first object so runs are reproducible.
    order = [v for v in objects[0].var_names if v in shared]

    if verbose:
        sizes = " / ".join(str(a.n_vars) for a in objects)
        print(f"gene-space intersection: {sizes} -> {len(order)} shared Ensembl IDs")

    return [adata[:, order].copy() for adata in objects]


# ---------------------------------------------------------------------------
# Back to symbols, once the merge is done
# ---------------------------------------------------------------------------

def to_canonical_symbols(
    adata,
    gene_space_dir: Path | None = None,
    copy: bool = False,
):
    """Switch `var_names` from Ensembl IDs to canonical (Ensembl 93) symbols.

    Call this ONCE, after the intersection/concat. At that point there is a single
    harmonized gene space, the mis-pairing risk is gone, and the ID's job is done —
    while every downstream consumer is symbol-native.

    The 9 symbols that still collide within the intersection (MATR3, RGS5, COG8, ...)
    are suffixed `SYMBOL__ENSGxxxxxxxxxxx`. That is deterministic and readable, unlike
    `var_names_make_unique()`, which assigns bare-vs-suffixed by row position and so
    reintroduces the exact ambiguity this module removes.

    Retains in `.var`: `ensembl_id`, `canonical_symbol`, `symbol_33538`,
    `symbol_33694`, `symbol_drift` (bool — did the symbol differ between builds).
    """
    if copy:
        adata = adata.copy()

    directory = gene_space_dir or config.GENE_SPACE_DIR
    path = directory / "gene_space_intersection.tsv.gz"
    if not path.exists():
        raise GeneSpaceError(f"Missing intersection table {path}.")

    table = pd.read_csv(path, sep="\t").set_index("ensembl_id")

    missing = [v for v in adata.var_names if v not in table.index]
    if missing:
        raise GeneSpaceError(
            f"{len(missing)} var_names are absent from the committed intersection "
            f"table (e.g. {missing[:3]}). Either the object was not intersected with "
            f"intersect_gene_space(), or the table is stale."
        )

    sub = table.loc[list(adata.var_names)]
    canonical = sub["canonical_symbol"].astype(str)

    duplicated = canonical.duplicated(keep=False)
    resolved = canonical.where(
        ~duplicated,
        canonical + "__" + pd.Series(adata.var_names, index=canonical.index).astype(str),
    )

    adata.var["ensembl_id"] = list(adata.var_names)
    adata.var["canonical_symbol"] = canonical.to_numpy()
    adata.var["symbol_33538"] = sub["symbol_33538"].to_numpy()
    adata.var["symbol_33694"] = sub["symbol_33694"].to_numpy()
    adata.var["symbol_drift"] = sub["symbol_drift"].astype(bool).to_numpy()
    # Index name is `symbol`, NOT `canonical_symbol`, and the difference is
    # load-bearing: `var["canonical_symbol"]` holds the BARE Ensembl-93 symbol while
    # the index holds the collision-resolved one, so for the 9 duplicated symbols the
    # two differ. AnnData refuses to write an index whose name matches a column with
    # different values, which made this an error that only appeared at `write_h5ad`
    # — long after the object looked correct in memory.
    adata.var_names = pd.Index(resolved.to_numpy(), name="symbol")

    if adata.var_names.duplicated().any():
        dupes = adata.var_names[adata.var_names.duplicated()].tolist()[:5]
        raise GeneSpaceError(f"Symbol disambiguation failed; still duplicated: {dupes}")

    return adata


# ---------------------------------------------------------------------------
# Assertions
# ---------------------------------------------------------------------------

def assert_required_genes(adata, groups: Iterable[str] | None = None) -> None:
    """Hard-fail with specific names if any required gene did not survive.

    These assertions caught the NSD2/WHSC1 drift that manual inspection had missed
    across two prior builds of this project. They stay, and they stay loud.

    A missing gene means "the join key regressed" or "check for a legacy symbol" —
    it does NOT mean the gene is biologically absent. All 65 required genes are
    confirmed present in the Ensembl-ID intersection.

    Works on an object keyed by either canonical symbols or Ensembl IDs (it reads
    `var["canonical_symbol"]` when present).
    """
    if "canonical_symbol" in adata.var:
        present = set(adata.var["canonical_symbol"].astype(str))
    else:
        present = set(map(str, adata.var_names))

    wanted = config.REQUIRED_GENES if groups is None else {
        g: config.REQUIRED_GENES[g] for g in groups
    }

    failures: list[str] = []
    for group, genes in wanted.items():
        missing = sorted(genes - present)
        if missing:
            hints = [
                f"{g} (legacy symbol {old!r} may be present instead)"
                for g in missing
                for old, new in config.LEGACY_SYMBOLS.items()
                if new == g
            ]
            failures.append(
                f"  {group}: {len(missing)} missing -> {missing}"
                + (f"\n    hint: {'; '.join(hints)}" if hints else "")
            )

    if failures:
        raise GeneSpaceError(
            "Required genes did not survive the gene-space intersection:\n"
            + "\n".join(failures)
            + "\n\nThis is a join-key or reference problem, not biology. Check that "
            "attach_ensembl_ids() ran on every object and that the intersection was "
            "on ensembl_id, not on symbols."
        )

    # Canary: if the ID join silently regressed to a symbol join, the drifted
    # canonical symbols are the first things to disappear.
    regressed = sorted(new for new in config.LEGACY_SYMBOLS.values() if new not in present)
    if regressed:
        raise GeneSpaceError(
            f"Canonical symbols {regressed} are absent while their legacy forms may "
            f"not be — the classic signature of a raw-symbol intersection. The join "
            f"must be on ensembl_id."
        )


# ---------------------------------------------------------------------------
# Regeneration (rarely needed — the map is committed)
# ---------------------------------------------------------------------------

def rebuild_gene_map_from_gtf(
    gtf_path: Path,
    build: int,
    deposited_symbols: Sequence[str],
) -> pd.DataFrame:
    """Reconstruct `deposited_symbol -> ensembl_id` from an Ensembl GTF, and verify.

    `deposited_symbols` is the literal `genes.tsv` column from any sample of that
    build (they are byte-identical within a build). The function raises unless the
    reconstruction reproduces it position-for-position — that assertion is the whole
    point, since it certifies the release and biotype filter were right.

    You should not normally need this: the map is committed at
    `resources/gene_space/`. Use it to re-derive from scratch, or to add a build.

        gtf = Path("Homo_sapiens.GRCh38.93.gtf.gz")   # config.BUILDS[33538]["gtf_url"]
        symbols = [l.strip() for l in open(".../genes.tsv")]
        frame = rebuild_gene_map_from_gtf(gtf, 33538, symbols)
    """
    rows: list[tuple[str, str]] = []
    seen: set[str] = set()

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            biotype = _GTF_GENE_BIOTYPE.search(attrs)
            if biotype is None or biotype.group(1) not in config.MKGTF_BIOTYPES:
                continue
            gene_id = _GTF_GENE_ID.search(attrs).group(1)
            gene_id = _VERSION_SUFFIX_RE.sub("", gene_id)
            if gene_id in seen:
                continue
            seen.add(gene_id)
            name = _GTF_GENE_NAME.search(attrs)
            rows.append((gene_id, name.group(1) if name else gene_id))

    if len(rows) != build:
        raise GeneSpaceError(
            f"GTF yielded {len(rows)} genes after the mkgtf biotype filter, expected "
            f"{build}. Wrong Ensembl release or wrong biotype whitelist — expected "
            f"release {config.BUILDS[build]['ensembl_release']}."
        )

    reconstructed = _make_unique([_seurat_sanitize(sym) for _, sym in rows])
    expected = list(deposited_symbols)
    if reconstructed != expected:
        mismatches = [
            (i, a, b) for i, (a, b) in enumerate(zip(reconstructed, expected)) if a != b
        ]
        raise GeneSpaceError(
            f"Reconstruction does not match the deposited column "
            f"({len(mismatches)} mismatches). First 5: {mismatches[:5]}"
        )

    return pd.DataFrame(
        {
            "deposited_symbol": expected,
            "ensembl_id": [gene_id for gene_id, _ in rows],
        },
        index=pd.RangeIndex(len(rows), name="row_index"),
    )
