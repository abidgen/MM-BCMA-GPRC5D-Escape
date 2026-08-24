"""
Loading GSE223060 from disk, and the sample metadata derivable from filenames.

WHY THIS MODULE EXISTS
----------------------
`scanpy.read_10x_mtx()` cannot read this deposit. The archive is old Cell Ranger
v2-style, uncompressed, one directory level deeper than the GSM folder, and the
matrix file is not called `matrix.mtx`:

    raw/samples/<GSM>_<sample>/<sample>/counts.mtx     (NOT matrix.mtx)
                                       /genes.tsv      (NOT features.tsv; ONE column)
                                       /barcodes.tsv   (uncompressed, no .gz)

`read_10x_mtx()` hardcodes the filenames it looks for, so pointing it at the right
directory fails to find anything. Everything here is built on explicit paths taken
from `raw/sample_manifest.csv` (produced by `scripts/03_build_manifest.py`, whose
column schema is the contract) instead of on any 10x-directory auto-detection.

`genes.tsv` is a single column of gene symbols — there is no `gene_id` column and no
`ENSG` string anywhere in the deposit. IDs are recovered separately and positionally
by `gene_space.attach_ensembl_ids`, which is why this module must return `var_names`
as the **deposited symbols in deposited order** and must never sort, dedupe or
otherwise touch the gene axis. `var_names_make_unique()` is banned for the same
reason — see `gene_space`.

ORIENTATION
-----------
The .mtx is genes x cells (Cell Ranger's orientation), e.g. `33538 1007 2604146`.
AnnData is cells x genes, so the matrix is transposed on read. Getting this backwards
produces an object that still looks plausible, which is why the dimensions are checked
against the two .tsv line counts rather than trusted.

PATIENT MAPPING — RESOLVED 2026-08-24 BY SUPPLEMENTARY TABLE S1
--------------------------------------------------------------
`patient_id` is now the S1-resolved mapping, and `obs["patient_id_source"]` reads
`"S1"` rather than `"naive"`. It is the naive rule (strip a trailing `_<digits>` only
when the stem is purely numeric) plus one alias — `MMY83942` is the same patient as
`83942` — and one sample, `25183`, that S1 does not list at all. Those two facts take
the deposit's 54 samples / 43 naive patients to the paper's 53 / 41, exactly; the
parser asserts both numbers. See `rebuild_clinical_metadata_from_s1`.

`naive_patient_id` is kept and still exported: it is what the assertion compares
against, and `patient_id_source` stays on every object so an aggregate built before
this change is still identifiable.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from anndata import AnnData

from . import config

__all__ = [
    "load_manifest",
    "load_sample_metadata",
    "load_clinical_metadata",
    "rebuild_sample_metadata_from_soft",
    "rebuild_clinical_metadata_from_s1",
    "s1_patient_id",
    "parse_sample_name",
    "naive_patient_id",
    "classify_sample",
    "read_sample",
    "read_samples",
    "SampleLoadError",
]


class SampleLoadError(RuntimeError):
    """Raised when a sample cannot be loaded as the pipeline requires.

    Always names the offending sample and file — a sample that loads into a
    plausible-looking but wrong object is the failure mode this module guards.
    """


# `GSM6939028_MMRF_1695` -> ("GSM6939028", "MMRF_1695")
_GSM_PREFIX_RE = re.compile(r"^(GSM\d+)_(.+)$")

#: Normal bone marrow controls. Confirmed by name convention; used by stage 07 as the
#: malignant-caller negative control (polyclonal marrow must yield no clone).
_NORMAL_BM_RE = re.compile(r"^BM\d+$")

#: `ND_083017`, `ND_090617`, `ND_170531`, `ND_170607` — normal donors. CONFIRMED
#: 2026-08-24 against the GEO SOFT file, which gives them `source_name` = "Donor BMMC,
#: aspirate, scRNAseq" and no `diagnosis` characteristic at all, while the other 54
#: samples read "Multiple myeloma (MM)". The suffixes are collection dates.
#:
#: This regex is now only a FALLBACK for when the committed metadata table is absent;
#: `load_manifest` prefers the GEO table, which is authoritative. It is kept because
#: the filename convention is the one thing that cannot go missing.
_NORMAL_DONOR_RE = re.compile(r"^ND_\d+$")

#: The naive patient rule. Strips a trailing `_<digits>` ONLY when the stem is purely
#: numeric: `27522_1` -> `27522`, but `MMRF_1695` stays whole (its stem is `MMRF`).
_NUMERIC_SUFFIX_RE = re.compile(r"^(\d+)_(\d+)$")


# ---------------------------------------------------------------------------
# Sample metadata derivable from names alone
# ---------------------------------------------------------------------------

def parse_sample_name(sample_id: str) -> tuple[str, str]:
    """Split a manifest `sample_id` into (gsm_id, sample_name).

    `GSM6939028_MMRF_1695` -> `("GSM6939028", "MMRF_1695")`. Splitting on the first
    underscore is not enough in general, so the GSM accession is matched explicitly.
    """
    match = _GSM_PREFIX_RE.match(sample_id)
    if match is None:
        raise SampleLoadError(
            f"sample_id {sample_id!r} does not look like <GSM accession>_<sample name>. "
            f"The manifest's sample_id column is the directory name under raw/samples/."
        )
    return match.group(1), match.group(2)


def naive_patient_id(sample_name: str) -> str:
    """Provisional patient ID — NOT the real mapping.

    Strips a trailing `_<digits>` only when the stem is purely numeric, so the
    six `27522_*` samples collapse to one patient while `MMRF_1695` is left alone.

    This yields 43 provisional patients from 53 disease samples against the paper's
    reported 41 / 53, so roughly two name collapses are still being missed
    (`83942` / `MMY83942` is the obvious candidate). Do NOT hand-fix those here —
    guessing a patient mapping silently changes the denominator of the headline
    metric. Resolve against Supplementary Table S1; see CLAUDE.md's S1 policy.

    The `_N` suffixes are also of *unknown meaning* (timepoint vs. fraction vs. sort
    vs. replicate) — the bulk/sc suffix misalignment is evidence against the naive
    timepoint reading. This function only assumes they belong to one patient.
    """
    match = _NUMERIC_SUFFIX_RE.match(sample_name)
    return match.group(1) if match else sample_name


def classify_sample(sample_name: str) -> tuple[str, bool]:
    """Return (sample_type, certain) for a sample name.

    Filename-only fallback, used when the committed GEO metadata table is missing.
    `load_manifest` overrides it from `resources/sample_metadata/`, which is
    authoritative and where `sample_type_certain` becomes True for every sample.
    """
    if _NORMAL_BM_RE.match(sample_name) or _NORMAL_DONOR_RE.match(sample_name):
        return "normal_bm", True
    return "myeloma", True


# ---------------------------------------------------------------------------
# GEO series metadata
# ---------------------------------------------------------------------------
#
# The SOFT files (raw/GSE22306{0,1}_family.soft.gz, ~8 KB each) carry per-sample
# facts that are NOT derivable from filenames and were not in the project's ground
# truth until 2026-08-24. Two of them are load-bearing:
#
#   cohort      MMRF / WU1 / WU2 / Donor
#   chemistry   10x 3' v2 (WashU cohort 1) vs v3.2 / v3.3 (everything else)
#
# Chemistry is confounded with cohort, and the resulting sensitivity difference is
# real but SMALLER THAN THE v2-vs-v3 FOLKLORE. Measured on this cohort's pre-QC cells
# (sample-level medians of genes detected per cell, 2026-08-24):
#
#     MMRF   v3.3   18 samples   1916 genes/cell
#     WU2    v3.2   13 samples   1210
#     Donor  v3.2    8 samples   1103
#     WU1    v2     23 samples   1023
#
#     v2 vs all-v3: 1023 vs 1408, a 1.38x ratio, Mann-Whitney p = 6.5e-05,
#     but the sample distributions OVERLAP (v2 max 1602 > v3 min 793).
#
# So the axis that actually separates is COHORT, not chemistry version per se — MMRF
# is ~1.9x the others and WU2/Donor sit close to WU1 despite being v3. Chemistry is
# one component of a cohort/site/protocol difference, not the whole of it.
#
# It still has to be modelled: the headline metric is a FRACTION OF ZEROS on a
# low-abundance transcript, so a 1.9x depth difference that tracks cohort will move
# frac_double_negative and read as biology. Carry `cohort` (and `chemistry`) as
# covariates in stage 08's depth regression and stage 10's null. Do NOT quote a
# "2-3x chemistry effect" — this cohort does not show one.
#
# `n_genes_ref` is NOT a usable proxy for it: the reference build cuts across cohorts
# (two WU1 samples on 33538, the four ND_* donors on 33694).
#
# raw/ is gitignored, so the parsed tables are committed under resources/ the way the
# gene-space map is, and `rebuild_sample_metadata_from_soft` regenerates them.

#: Prep protocol per cohort, read from !Sample_extract_protocol_ch1. WashU cohort 1
#: was loaded straight after thaw; every other cohort went through Miltenyi dead-cell
#: removal first, which is a second, smaller batch axis on top of chemistry.
COHORT_PROTOCOL: dict[str, dict[str, str]] = {
    "MMRF":  {"chemistry": "10x 3' v3.3", "dead_cell_removal": "yes"},
    "WU1":   {"chemistry": "10x 3' v2",   "dead_cell_removal": "no"},
    "WU2":   {"chemistry": "10x 3' v3.2", "dead_cell_removal": "yes"},
    "Donor": {"chemistry": "10x 3' v3.2", "dead_cell_removal": "yes"},
}

_SOURCE_COHORT = (
    ("MMRF cohort", "MMRF"),
    ("WU cohort 1", "WU1"),
    ("WU cohort 2", "WU2"),
    ("Donor", "Donor"),
)


def _soft_records(soft_path: Path) -> list[dict[str, str]]:
    """Parse a GEO `*_family.soft.gz` into one flat dict per ^SAMPLE block."""
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    opener = gzip.open if str(soft_path).endswith(".gz") else open
    with opener(soft_path, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                current = {"gsm_id": line.split(" = ", 1)[1].strip()}
                records.append(current)
            elif current is None or not line.startswith("!Sample_"):
                continue
            elif " = " in line:
                key, value = line[len("!Sample_"):].split(" = ", 1)
                if key == "characteristics_ch1" and ": " in value:
                    sub, value = value.split(": ", 1)
                    key = f"char_{sub.strip().replace(' ', '_')}"
                # Repeated keys (data_processing, description) are joined, not lost.
                current[key] = f"{current[key]} || {value}" if key in current else value
    return records


def rebuild_sample_metadata_from_soft(
    soft_path: Path,
    assay: str,
    out_path: Path | None = None,
) -> pd.DataFrame:
    """Parse a GEO SOFT file into the committed per-sample metadata table.

    `assay` is `"scrna"` (GSE223060) or `"bulk"` (GSE223061). Emits `gsm_id`,
    `sample_name`, `cohort`, `diagnosis`, `sample_type`, `source_name`, plus
    `chemistry`/`dead_cell_removal` for scRNA and `prep` for bulk.

    The committed tables live at `resources/sample_metadata/`; this only needs
    re-running if GEO revises the deposit. Like the gene-space reconstruction, it
    asserts its own expectations (sample count, that every sample resolves to a
    known cohort) so a changed deposit fails loudly rather than merging silently.
    """
    records = _soft_records(soft_path)
    expected = {"scrna": 62, "bulk": 31}[assay]
    if len(records) != expected:
        raise SampleLoadError(
            f"{soft_path} holds {len(records)} samples, expected {expected} for "
            f"assay={assay!r}. The deposit changed — re-verify before using it."
        )

    rows: list[dict[str, object]] = []
    for record in records:
        source = record.get("source_name_ch1", "")
        cohort = next((c for token, c in _SOURCE_COHORT if token in source), None)
        if cohort is None:
            raise SampleLoadError(
                f"{record['gsm_id']}: source_name {source!r} matches no known cohort."
            )
        name = record.get("title", "").replace("bulk_RNA_", "")
        # Absent `diagnosis` IS the control marker: the 8 donor samples carry no
        # diagnosis characteristic at all, the other 54 read "Multiple myeloma (MM)".
        diagnosis = record.get("char_diagnosis", "")
        row: dict[str, object] = {
            "gsm_id": record["gsm_id"],
            "sample_name": name,
            "cohort": cohort,
            "diagnosis": diagnosis or "none",
            "sample_type": "normal_bm" if cohort == "Donor" else "myeloma",
            "source_name": source,
        }
        if assay == "scrna":
            row.update(COHORT_PROTOCOL[cohort])
        else:
            row["prep"] = (
                "CD138+ sorted" if "CD138+ sorted" in source else "unsorted BMMC"
            )
        rows.append(row)

    frame = pd.DataFrame(rows)

    disagree = frame.loc[
        (frame["sample_type"] == "myeloma") != frame["diagnosis"].str.contains("myeloma"),
        "sample_name",
    ].tolist()
    if disagree:
        raise SampleLoadError(
            f"cohort and diagnosis disagree for {disagree}. Donor samples must carry "
            f"no diagnosis and MM samples must carry one."
        )

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out_path, sep="\t", index=False)
    return frame


def load_sample_metadata(assay: str = "scrna", directory: Path | None = None) -> pd.DataFrame:
    """Load the committed GEO metadata table for `"scrna"` or `"bulk"`."""
    directory = directory or config.SAMPLE_METADATA_DIR
    path = directory / f"{assay}_samples.tsv"
    if not path.exists():
        raise SampleLoadError(
            f"Missing {path}. Regenerate with rebuild_sample_metadata_from_soft("
            f"raw/GSE223060_family.soft.gz, {assay!r}, {path})."
        )
    return pd.read_csv(path, sep="\t")


# ---------------------------------------------------------------------------
# Supplementary Table S1 — the clinical metadata, and the real patient mapping
# ---------------------------------------------------------------------------

#: S1 writes cohorts out longhand; the rest of the project uses the GEO vocabulary.
_S1_COHORT = {
    "WASHU Cohort 1": "WU1",
    "WASHU Cohort 2": "WU2",
    "MMRF Immune Atlas Pilot Study Cohort": "MMRF",
}

#: The disease-stage sheet spells a per-patient serial course: `Primary`, `SMM`, and
#: numbered `Relapse-N`/`Remission-N`. Collapsed to a coarse phase because n is small
#: and the numbered levels are per patient, not comparable across them.
_DISEASE_PHASE = {
    "SMM": "smoldering",
    "Primary": "newly_diagnosed",
    "Relapse": "relapsed",
    "Remission": "remission",
}


def s1_patient_id(sample_name: str) -> str:
    """The canonical patient id for a sample name, S1 aliases applied.

    `naive_patient_id` plus `config.PATIENT_ALIASES`. The only alias is
    `MMY83942` -> `83942`; see the constant for the evidence and the arithmetic.
    """
    naive = naive_patient_id(sample_name)
    return config.PATIENT_ALIASES.get(naive, naive)


def _s1_sheet(xlsx_path: Path, sheet: str) -> pd.DataFrame:
    try:
        return pd.read_excel(xlsx_path, sheet_name=sheet, header=None)
    except Exception as error:  # openpyxl raises several unrelated types
        raise SampleLoadError(
            f"Could not read sheet {sheet!r} from {xlsx_path}: {error}"
        ) from error


def rebuild_clinical_metadata_from_s1(
    xlsx_path: Path | None = None,
    directory: Path | None = None,
    *,
    write: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse Supplementary Table S1 into the two committed clinical tables.

    Returns `(patients, sample_stages)` and, unless `write=False`, commits them to
    `resources/sample_metadata/patients_clinical.tsv` and `sample_disease_stage.tsv`
    — same pattern as the GEO tables, so `raw/` staying gitignored costs nothing.

    WHAT S1 SETTLES
    ---------------
    Two things that had been open since the project started:

    1. **The patient mapping.** The naive rule gives 43 patients over 54 deposited
       myeloma samples; the paper reports 41 over 53. Both gaps close here and
       exactly: `25183` is deposited but appears in no supplementary table
       (`config.SAMPLES_WITHOUT_CLINICAL`), and `83942`/`MMY83942` are one patient
       sampled under both WashU protocols (`config.PATIENT_ALIASES`). 54 - 1 = 53
       samples; 43 - 1 - 1 = 41 patients. The function asserts both numbers.

    2. **What the `_N` suffixes mean.** Sheet 2 reads `27522_1` Primary, `_2`
       Remission-1, `_3` Relapse-1, `_4` Relapse-2, `_5` Remission-2, `_6`
       Relapse-3. They are serial disease-course timepoints, which is what the
       2026-08-24 bulk/scRNA suffix argument had inferred and this confirms
       outright. The longitudinal arm in CLAUDE.md is real, not speculative.
       It also explains the lone non-`_1` samples: `37692_2` and `57075_3` are later
       timepoints whose earlier draws were simply not deposited.

    COVERAGE, STATED PLAINLY
    ------------------------
    Sheet 2 covers **WashU cohort 1 only**. MMRF and WashU cohort 2 samples get
    `disease_stage = NA` — S1 gives them treatment and time-to-progression but no
    stage label. Do not impute one; `newly_diagnosed` is a guess for those cohorts,
    not a datum.
    """
    xlsx_path = Path(xlsx_path) if xlsx_path is not None else config.S1_XLSX
    if not xlsx_path.exists():
        raise SampleLoadError(
            f"Supplementary Table S1 not found at {xlsx_path}. It is not in git "
            f"(raw/ is gitignored); download it from the Cancer Research page for "
            f"CAN-22-1769 as can-22-1769_table_s1_suppst1.xlsx."
        )

    # --- sheet 1: per-patient clinical characteristics -----------------------
    raw = _s1_sheet(xlsx_path, "summary")
    # Row 0 is the table title, row 1 the header, rows 2+ the data. Column 0 names the
    # cohort on its first patient only, so it is forward-filled.
    body = raw.iloc[2:, :8].copy()
    body.columns = [
        "s1_cohort", "patient_id", "age", "sex", "race",
        "iss_stage", "treatment", "ttpd_months",
    ]
    body["s1_cohort"] = body["s1_cohort"].ffill()
    body = body.loc[body["patient_id"].notna()].reset_index(drop=True)

    unknown_cohorts = sorted(set(body["s1_cohort"]) - set(_S1_COHORT))
    if unknown_cohorts:
        raise SampleLoadError(
            f"S1 sheet 'summary' names cohort(s) this parser does not know: "
            f"{unknown_cohorts}. Extend _S1_COHORT rather than guessing."
        )
    body["cohort"] = body["s1_cohort"].map(_S1_COHORT)
    body["patient_id"] = body["patient_id"].astype(str).str.strip()

    # ISS is 1/2/3 or the literal "UNK"; keep it as a string so "UNK" is not coerced
    # to NaN and silently read as missing-at-random.
    body["iss_stage"] = body["iss_stage"].astype(str).str.strip().replace(
        {"nan": "UNK", "": "UNK"}
    )
    body["ttpd_months"] = pd.to_numeric(body["ttpd_months"], errors="coerce")
    body["age"] = pd.to_numeric(body["age"], errors="coerce").astype("Int64")
    body["treatment"] = body["treatment"].fillna("Unknown").astype(str).str.strip()

    # Fold the alias BEFORE the count assertion — that collapse is the point.
    body["s1_listed_as"] = body["patient_id"]
    body["patient_id"] = body["patient_id"].map(
        lambda pid: config.PATIENT_ALIASES.get(pid, pid)
    )
    aliased = sorted(set(body["s1_listed_as"]) & set(config.PATIENT_ALIASES))
    # Two rows now share a patient_id. Keep the first (the cohort-1 row) and record
    # that the patient spans two cohorts rather than dropping the fact.
    spans = body.groupby("patient_id")["cohort"].apply(lambda c: "+".join(sorted(set(c))))
    patients = body.drop_duplicates("patient_id", keep="first").set_index("patient_id")
    patients["cohorts_sampled"] = spans
    patients = patients.reset_index()

    # --- sheet 2: per-sample disease stage (WashU cohort 1 only) -------------
    stage_raw = _s1_sheet(xlsx_path, "WashU cohort1 stage info")
    stages = stage_raw.iloc[1:, :2].copy()
    stages.columns = ["sample_name", "disease_stage"]
    stages = stages.loc[stages["sample_name"].notna()].reset_index(drop=True)
    stages["sample_name"] = stages["sample_name"].astype(str).str.strip()
    stages["disease_stage"] = stages["disease_stage"].astype(str).str.strip()
    # `Relapse-2` -> `relapsed`; `Primary` -> `newly_diagnosed`; `SMM` -> `smoldering`.
    stages["disease_phase"] = (
        stages["disease_stage"].str.split("-").str[0].map(_DISEASE_PHASE)
    )
    unmapped = sorted(set(stages.loc[stages["disease_phase"].isna(), "disease_stage"]))
    if unmapped:
        raise SampleLoadError(
            f"S1 disease stages this parser cannot phase: {unmapped}. "
            f"Extend _DISEASE_PHASE rather than dropping them."
        )
    stages["patient_id"] = stages["sample_name"].map(s1_patient_id)
    # Timepoint index straight off the deposited suffix; NA for unsuffixed samples.
    stages["timepoint"] = (
        stages["sample_name"].str.extract(r"_(\d+)$")[0].astype("Int64")
    )

    _assert_s1_reproduces_the_paper(patients, aliased)

    if write:
        directory = directory or config.SAMPLE_METADATA_DIR
        directory.mkdir(parents=True, exist_ok=True)
        patients.to_csv(directory / "patients_clinical.tsv", sep="\t", index=False)
        stages.to_csv(directory / "sample_disease_stage.tsv", sep="\t", index=False)

    return patients, stages


def _assert_s1_reproduces_the_paper(
    patients: pd.DataFrame, aliased: list[str]
) -> None:
    """Check the S1 parse against the deposit and the paper's stated counts.

    This is the same self-certifying discipline as the gene-space reconstruction: the
    mapping is only trustworthy because getting it wrong changes a number that is
    checked here. It runs off the committed GEO metadata table, so it needs no `raw/`
    matrices.
    """
    if not aliased:
        raise SampleLoadError(
            f"Expected S1 to list the aliased patient(s) "
            f"{sorted(config.PATIENT_ALIASES)}, but none appeared. Either S1 was "
            f"revised or config.PATIENT_ALIASES is stale — do not proceed on a "
            f"patient mapping that no longer matches its evidence."
        )

    meta = load_sample_metadata("scrna")
    myeloma = meta.loc[meta["sample_type"] == "myeloma", "sample_name"]
    if len(myeloma) != config.N_MYELOMA_SAMPLES_DEPOSITED:
        raise SampleLoadError(
            f"GEO metadata holds {len(myeloma)} myeloma samples, expected "
            f"{config.N_MYELOMA_SAMPLES_DEPOSITED}."
        )

    in_paper = myeloma[~myeloma.isin(config.SAMPLES_WITHOUT_CLINICAL)]
    if len(in_paper) != config.N_MYELOMA_SAMPLES_IN_PAPER:
        raise SampleLoadError(
            f"{len(in_paper)} myeloma samples carry clinical data, expected "
            f"{config.N_MYELOMA_SAMPLES_IN_PAPER} (the paper's 53). "
            f"config.SAMPLES_WITHOUT_CLINICAL may be stale."
        )

    resolved = {s1_patient_id(name) for name in in_paper}
    if len(resolved) != config.N_PATIENTS_IN_PAPER:
        raise SampleLoadError(
            f"The S1 mapping yields {len(resolved)} patients over the deposited "
            f"myeloma samples, expected {config.N_PATIENTS_IN_PAPER}. The mapping "
            f"sets the denominator of frac_double_negative — fix it, do not proceed."
        )

    # Every deposited patient must have a clinical row, or the join below is partial.
    missing = sorted(resolved - set(patients["patient_id"]))
    if missing:
        raise SampleLoadError(
            f"{len(missing)} deposited patient(s) have no S1 clinical row: {missing}. "
            f"Add them to config.SAMPLES_WITHOUT_CLINICAL only with evidence."
        )


def load_clinical_metadata(
    directory: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the committed `(patients, sample_stages)` clinical tables from S1."""
    directory = directory or config.SAMPLE_METADATA_DIR
    paths = {
        name: directory / f"{name}.tsv"
        for name in ("patients_clinical", "sample_disease_stage")
    }
    absent = [str(p) for p in paths.values() if not p.exists()]
    if absent:
        raise SampleLoadError(
            f"Missing clinical table(s) {absent}. Regenerate with "
            f"rebuild_clinical_metadata_from_s1()."
        )
    return tuple(pd.read_csv(p, sep="\t") for p in paths.values())  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# The manifest
# ---------------------------------------------------------------------------

def load_manifest(
    path: Path | None = None,
    *,
    drop_excluded: bool = True,
    require_complete: bool = True,
    with_metadata: bool = True,
    with_clinical: bool = True,
) -> pd.DataFrame:
    """Load `raw/sample_manifest.csv` and add the metadata derivable from names.

    The manifest stores **repo-root-relative** paths (it is committed and read on
    other machines), so they are resolved to absolute here and their existence is
    checked. Columns added on top of the script's schema:

        gsm_id, sample_name, patient_id, patient_id_naive, patient_id_source,
        sample_type, sample_type_certain, excluded

    from the committed GEO metadata table (`with_metadata`, the default):

        cohort, chemistry, dead_cell_removal, diagnosis

    and from the committed Supplementary Table S1 tables (`with_clinical`, the
    default):

        age, sex, race, iss_stage, treatment, ttpd_months,
        disease_stage, disease_phase, timepoint,
        clinical_source, in_paper_cohort

    `disease_stage`/`disease_phase`/`timepoint` are **WashU cohort 1 only** — S1
    gives no stage label for MMRF or WashU cohort 2, and one is not imputed. Where
    present they are serial timepoints (`27522_1` Primary -> `_6` Relapse-3), which
    is what makes the longitudinal arm real rather than assumed.

    `cohort` is the one to watch. WashU cohort 1 ran 10x 3' v2 and everything else ran
    v3.2/v3.3, but the measured gap in genes detected per cell follows cohort rather
    than chemistry version: MMRF 1916, WU2 1210, Donor 1103, WU1 1023 (sample-level
    medians, pre-QC). v2-vs-v3 is 1.38x with overlapping distributions, not the 2-3x
    the chemistry difference might suggest. Since the headline metric is a fraction of
    zeros, a 1.9x depth spread across cohorts must be carried as a covariate.

    `drop_excluded` filters `config.EXCLUDED_SAMPLES`, which is currently **empty** —
    `56203_1` was excluded until 2026-08-24 on a misdiagnosis and is now repaired on
    read instead (see `config.TRUNCATED_GENE_FILES`). The parameter stays for the next
    sample that genuinely has to go.

    `require_complete` hard-fails on any row not classified `triplet-ok` by the
    manifest script, rather than discovering the problem mid-load.

    NOTE on counts: 62 rows, of which **8 are normal-BM controls** (`BM2/4/5/6` plus
    the four `ND_*`) and **54 are myeloma** — confirmed against GEO, not inferred.
    Under `s1_patient_id` the 54 resolve to **41 patients over 53 in-cohort samples**,
    matching the paper: `25183` carries no S1 entry (`in_paper_cohort == False`, but
    it is still loaded) and `MMY83942` folds into `83942`. The naive rule's 43 is
    retained as `patient_id_naive` for comparison. An earlier figure of "47 patients /
    57 samples" counted the four donors as disease and is superseded.
    """
    manifest_path = Path(path) if path is not None else config.MANIFEST_CSV
    if not manifest_path.exists():
        raise SampleLoadError(
            f"No manifest at {manifest_path}. Build it with "
            f"`python scripts/03_build_manifest.py raw/samples` or "
            f"`notebooks/03_build_manifest.py`."
        )

    frame = pd.read_csv(manifest_path)

    if require_complete:
        bad = frame.loc[frame["format"] != "triplet-ok", "sample_id"].tolist()
        if bad:
            raise SampleLoadError(
                f"{len(bad)} sample(s) are not 'triplet-ok' in {manifest_path}: {bad}. "
                f"This loader reads the counts.mtx/genes.tsv/barcodes.tsv triplet only; "
                f"re-run the manifest builder or inspect those directories."
            )

    parsed = frame["sample_id"].map(parse_sample_name)
    frame["gsm_id"] = [gsm for gsm, _ in parsed]
    frame["sample_name"] = [name for _, name in parsed]
    frame["patient_id"] = frame["sample_name"].map(s1_patient_id)
    frame["patient_id_source"] = "S1"
    frame["patient_id_naive"] = frame["sample_name"].map(naive_patient_id)

    classified = frame["sample_name"].map(classify_sample)
    frame["sample_type"] = [t for t, _ in classified]
    frame["sample_type_certain"] = [c for _, c in classified]

    # The GEO metadata is authoritative where it exists; the filename rules above are
    # the fallback. Joining on sample_name rather than gsm_id keeps this working if a
    # sample is ever re-accessioned.
    if with_metadata:
        meta = load_sample_metadata("scrna").set_index("sample_name")
        unknown = sorted(set(frame["sample_name"]) - set(meta.index))
        if unknown:
            raise SampleLoadError(
                f"{len(unknown)} sample(s) are on disk but absent from the committed "
                f"GEO metadata table: {unknown[:5]}. Regenerate it with "
                f"rebuild_sample_metadata_from_soft()."
            )
        joined = meta.loc[frame["sample_name"]]
        for column in ("cohort", "chemistry", "dead_cell_removal", "diagnosis"):
            frame[column] = joined[column].to_numpy()
        # GEO overrides the filename guess, and settles it.
        frame["sample_type"] = joined["sample_type"].to_numpy()
        frame["sample_type_certain"] = True

    # Clinical metadata from S1. Left-joined on the resolved patient id, so the 8
    # donors and `25183` come through with NA rather than being dropped — the donors
    # are stage 07's negative control and must survive every join in this module.
    if with_clinical:
        patients, stages = load_clinical_metadata()
        clinical = patients.set_index("patient_id")
        for column in ("age", "sex", "race", "iss_stage", "treatment", "ttpd_months"):
            frame[column] = frame["patient_id"].map(clinical[column])
        stage_by_sample = stages.set_index("sample_name")
        frame["disease_stage"] = frame["sample_name"].map(
            stage_by_sample["disease_stage"]
        )
        frame["disease_phase"] = frame["sample_name"].map(
            stage_by_sample["disease_phase"]
        )
        frame["timepoint"] = frame["sample_name"].map(stage_by_sample["timepoint"])
        frame["clinical_source"] = np.where(
            frame["patient_id"].isin(clinical.index), "S1", "none"
        )
        frame["in_paper_cohort"] = (frame["sample_type"] == "myeloma") & ~frame[
            "sample_name"
        ].isin(config.SAMPLES_WITHOUT_CLINICAL)

    frame["excluded"] = frame["sample_name"].isin(config.EXCLUDED_SAMPLES)

    # The manifest holds repo-root-relative paths; resolve and verify them here so a
    # missing file surfaces now rather than 40 samples into a loop.
    root = config.REPO_ROOT
    for column in ("matrix_path", "barcodes_path", "genefeat_path"):
        frame[column] = [
            "" if not value else str((root / value).resolve())
            for value in frame[column].fillna("")
        ]

    if drop_excluded:
        frame = frame.loc[~frame["excluded"]].reset_index(drop=True)

    missing = [
        (row.sample_id, column)
        for row in frame.itertuples()
        for column in ("matrix_path", "barcodes_path", "genefeat_path")
        if not Path(getattr(row, column)).exists()
    ]
    if missing:
        raise SampleLoadError(
            f"{len(missing)} manifest path(s) do not exist on disk, e.g. {missing[:3]}. "
            f"Was raw/ extracted on this machine? See notebooks/01_download_data.py."
        )

    return frame


# ---------------------------------------------------------------------------
# Reading one sample
# ---------------------------------------------------------------------------

def _read_single_column(path: Path, what: str, sample_id: str) -> list[str]:
    """Read a one-column .tsv into a list of stripped strings.

    `genes.tsv` here is genuinely one column (no Ensembl ID), and `barcodes.tsv` is
    too. A file that turns out to have more columns is a different deposit than the
    one this pipeline was verified against, so it raises instead of silently taking
    the first field.
    """
    values: list[str] = []
    with open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            value = line.rstrip("\n").rstrip("\r")
            if "\t" in value:
                raise SampleLoadError(
                    f"{sample_id}: {what} file {path} has more than one column at line "
                    f"{lineno} ({value!r}). This deposit's {what} file is single-column; "
                    f"a multi-column file means a different reference format."
                )
            values.append(value)
    if not values:
        raise SampleLoadError(f"{sample_id}: {what} file {path} is empty.")
    return values


def _repair_truncated_genes(sample_name: str, truncated: list[str]) -> list[str]:
    """Replace a truncated `genes.tsv` column with the canonical one for its build.

    `56203_1`'s gene file stops mid-symbol at row 22185 (`KBTBD`, where the reference
    has `KBTBD7`) while its matrix declares the full 33694 rows — a failed write, not
    a different reference. The repair substitutes the canonical column from the
    committed, position-verified gene map.

    What makes this a repair rather than a guess is the prefix assertion: every row
    the deposit *did* write must match the canonical column exactly, and the final
    partial row must be a prefix of the symbol it was cut out of. If either fails,
    the file is damaged in some other way and this raises instead of substituting.
    """
    from . import gene_space  # local: gene_space imports config, io does not import it

    spec = config.TRUNCATED_GENE_FILES[sample_name]
    build = spec["build"]
    if len(truncated) != spec["deposited_rows"]:
        raise SampleLoadError(
            f"{sample_name}: genes.tsv has {len(truncated)} rows, but the recorded "
            f"truncation is at {spec['deposited_rows']}. The file on disk is not the "
            f"one this repair was verified against — re-verify before loading it."
        )

    canonical = list(gene_space.load_gene_map(build)["deposited_symbol"])
    intact, partial = truncated[:-1], truncated[-1]
    if intact != canonical[: len(intact)]:
        first = next(
            i for i, (a, b) in enumerate(zip(intact, canonical)) if a != b
        )
        raise SampleLoadError(
            f"{sample_name}: genes.tsv is not a prefix of the {build}-gene reference "
            f"(diverges at row {first + 1}: {intact[first]!r} vs "
            f"{canonical[first]!r}). It is damaged in some way other than truncation; "
            f"do not substitute."
        )
    if not canonical[len(intact)].startswith(partial):
        raise SampleLoadError(
            f"{sample_name}: the final written row {partial!r} is not a prefix of "
            f"{canonical[len(intact)]!r}, so the file did not simply stop mid-symbol."
        )

    print(
        f"  {sample_name}: repaired truncated genes.tsv "
        f"({len(truncated)} written rows -> {len(canonical)}, "
        f"prefix verified against the {build}-gene reference)"
    )
    return canonical


def read_sample(
    sample: str | pd.Series,
    manifest: pd.DataFrame | None = None,
    *,
    dtype: str = "float32",
) -> AnnData:
    """Read one sample's counts.mtx/genes.tsv/barcodes.tsv triplet into an AnnData.

    `sample` is either a manifest `sample_id` (`GSM6939028_MMRF_1695`), a bare
    `sample_name` (`MMRF_1695`), or a manifest row. Passing `manifest` avoids
    re-reading the CSV in a loop.

    Returns cells x genes with:
        X            raw integer counts, CSR, as `dtype` (float32 by default — the
                     values stay integral, but int64 would triple the memory and
                     scanpy normalizes to float on the first pp step anyway)
        var_names    the DEPOSITED symbols in DEPOSITED ORDER, untouched. This is a
                     precondition of gene_space.attach_ensembl_ids, which verifies
                     them position-for-position.
        obs_names    `<sample_name>_<barcode>` — barcodes repeat across samples, so
                     they must be disambiguated before any concat. The original is
                     kept in obs["barcode"].
        obs          sample_id, gsm_id, sample_name, patient_id, patient_id_source,
                     sample_type, sample_type_certain, n_genes_ref, genes_repaired,
                     and (from the GEO table) cohort, chemistry, dead_cell_removal,
                     diagnosis
        uns          the source paths and the reference row count

    `n_genes_ref` is carried in `.obs` (not only `.uns`) because it survives
    concatenation and stage 05 needs it as a Harmony covariate.
    """
    if manifest is None:
        manifest = load_manifest()

    if isinstance(sample, pd.Series):
        row = sample
    else:
        hits = manifest.loc[
            (manifest["sample_id"] == sample) | (manifest["sample_name"] == sample)
        ]
        if len(hits) == 0:
            if sample in config.EXCLUDED_SAMPLES:
                raise SampleLoadError(
                    f"{sample!r} is on config.EXCLUDED_SAMPLES and was dropped from "
                    f"the manifest. Pass drop_excluded=False to inspect it."
                )
            raise SampleLoadError(f"{sample!r} matches no manifest row.")
        if len(hits) > 1:
            raise SampleLoadError(
                f"{sample!r} matches {len(hits)} manifest rows "
                f"({hits['sample_id'].tolist()}); pass the full sample_id."
            )
        row = hits.iloc[0]

    sample_id = str(row["sample_id"])
    if bool(row.get("excluded", False)):
        raise SampleLoadError(
            f"{sample_id} is excluded (see config.EXCLUDED_SAMPLES) and must not be "
            f"loaded for analysis."
        )

    matrix_path = Path(row["matrix_path"])
    genes_path = Path(row["genefeat_path"])
    barcodes_path = Path(row["barcodes_path"])

    symbols = _read_single_column(genes_path, "genes", sample_id)
    barcodes = _read_single_column(barcodes_path, "barcodes", sample_id)

    repaired = False
    if str(row['sample_name']) in config.TRUNCATED_GENE_FILES:
        symbols = _repair_truncated_genes(str(row['sample_name']), symbols)
        repaired = True

    # spmatrix=True is explicit, not decorative: scipy 1.20 flips the default to
    # sparse *arrays*, and the scanpy/anndata stack this feeds is still spmatrix-native.
    # Pinning it here keeps the loader's return type stable across that change.
    matrix = scipy.io.mmread(matrix_path, spmatrix=True)
    if matrix.shape != (len(symbols), len(barcodes)):
        raise SampleLoadError(
            f"{sample_id}: {matrix_path.name} is {matrix.shape[0]}x{matrix.shape[1]} but "
            f"genes.tsv has {len(symbols)} rows and barcodes.tsv has {len(barcodes)}. "
            f"Cell Ranger writes genes x cells; a swap here would transpose the whole "
            f"object into something that still looks plausible."
        )

    # genes x cells -> cells x genes. CSR after the transpose, not before.
    counts = sp.csr_matrix(matrix.T, dtype=np.dtype(dtype))

    # The deposited column already went through R's make.unique, so duplicates mean
    # the file is not what gene_space's positional reconstruction was verified against.
    if len(set(symbols)) != len(symbols):
        duplicated = sorted({s for s in symbols if symbols.count(s) > 1})[:5]
        raise SampleLoadError(
            f"{sample_id}: genes.tsv has duplicate symbols (e.g. {duplicated}). The "
            f"deposited column is make.unique'd and must be unique. Do NOT call "
            f"var_names_make_unique() to paper over this — see gene_space."
        )
    if len(set(barcodes)) != len(barcodes):
        raise SampleLoadError(f"{sample_id}: barcodes.tsv has duplicate barcodes.")

    sample_name = str(row["sample_name"])
    obs = pd.DataFrame(
        {
            "barcode": barcodes,
            "sample_id": sample_id,
            "gsm_id": str(row["gsm_id"]),
            "sample_name": sample_name,
            "patient_id": str(row["patient_id"]),
            "patient_id_source": str(row["patient_id_source"]),
            "sample_type": str(row["sample_type"]),
            "sample_type_certain": bool(row["sample_type_certain"]),
            **{
                column: str(row[column])
                for column in ("cohort", "chemistry", "dead_cell_removal", "diagnosis")
                if column in row
            },
            "n_genes_ref": len(symbols),
            "genes_repaired": repaired,
        },
        index=pd.Index([f"{sample_name}_{bc}" for bc in barcodes], name="cell_id"),
    )
    for column in (
        "sample_id", "gsm_id", "sample_name", "patient_id", "patient_id_source",
        "sample_type", "n_genes_ref", "cohort", "chemistry", "dead_cell_removal",
        "diagnosis",
    ):
        if column in obs:
            obs[column] = obs[column].astype("category")

    adata = AnnData(
        X=counts,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(symbols, name="deposited_symbol")),
    )
    adata.uns["source"] = {
        "matrix": str(matrix_path),
        "genes": str(genes_path),
        "barcodes": str(barcodes_path),
        "n_genes_ref": len(symbols),
    }
    adata.uns["genes_repaired"] = repaired
    return adata


def read_samples(
    samples: Iterable[str] | None = None,
    manifest: pd.DataFrame | None = None,
    *,
    verbose: bool = True,
    **kwargs,
) -> Iterator[AnnData]:
    """Yield one AnnData per sample, in manifest order.

    A generator on purpose: the full cohort is 204,040 pre-QC cells across 62 samples
    and stage 04 checkpoints each sample's post-QC object individually, so nothing
    needs all 62 raw matrices resident at once.
    """
    if manifest is None:
        manifest = load_manifest()
    if samples is not None:
        wanted = set(samples)
        manifest = manifest.loc[
            manifest["sample_id"].isin(wanted) | manifest["sample_name"].isin(wanted)
        ]
        found = set(manifest["sample_id"]) | set(manifest["sample_name"])
        unknown = sorted(wanted - found)
        if unknown:
            raise SampleLoadError(f"Unknown sample(s): {unknown}")

    total = len(manifest)
    for i, (_, row) in enumerate(manifest.iterrows(), start=1):
        adata = read_sample(row, manifest=manifest, **kwargs)
        if verbose:
            print(
                f"[{i:>2}/{total}] {row['sample_id']:<26} "
                f"{adata.n_obs:>5} cells x {adata.n_vars} genes  "
                f"({row['sample_type']})"
            )
        yield adata
