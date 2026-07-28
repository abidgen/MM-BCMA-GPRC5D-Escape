# What We're Doing, and Why — Plain-Language Overview (Revised)

**Revision note:** The original version of this doc assumed raw FASTQ data was available
for GSE223060/GSE223061 and planned to process it through NCI's SINCLAIR pipeline on
Biowulf. That assumption was wrong — confirmed directly against NCBI: this BioProject
(PRJNA924769) has no SRA accession and no raw sequencing reads deposited anywhere
publicly. What GEO hosts instead is Cell Ranger's *output* — per-sample count matrices —
bundled as `GSE223060_RAW.tar` (~970 MB, 62 samples). This doc reflects that reality.
Everything from Stage B onward is unchanged from the original plan.

---

## 1. The business problem this maps to

Legend Biotech's core product, CARVYKTI, is a CAR-T therapy that targets **BCMA**, a
protein on multiple myeloma cells. It works well, but a meaningful share of patients
relapse — and one of the known reasons is **antigen escape**: the myeloma cells that
survive treatment are the ones that stopped expressing BCMA. The field's response has
been to develop a second target, **GPRC5D**, either as a fallback after BCMA fails or
paired with BCMA in a dual-target construct.

The open question this analysis speaks to: **if you targeted both BCMA and GPRC5D at
once, how many patients would still have some malignant cells slip through because
those cells don't express *either* antigen?** That's not a question the antigen-escape
literature answers directly — most papers ask "did this patient's antigen disappear
after treatment," which is a *before/after* question. We're asking a *right-now,
pre-treatment* question: how much of the risk is already baked into the tumor's
existing diversity, before treatment even starts.

---

## 2. The dataset — what's actually available

**GSE223060 / GSE223061** — single-cell RNA-seq and matched bulk RNA-seq from bone
marrow samples across multiple myeloma patients (three cohorts: the MMRF Immune Atlas
Pilot study and two internal WashU cohorts), generated specifically to discover new
therapeutic surface targets in myeloma.

**What's public:** processed, per-sample Cell Ranger output — 62 sample bundles,
~970 MB total — not raw FASTQ. There is no SRA accession for this data; it was never
deposited as raw reads. This is a real constraint, not a technical hiccup: it means the
alignment step (raw reads → per-cell gene counts) was already done by the original
authors, and you're starting one step downstream of that.

scRNA-seq still means what it always meant: instead of one averaged gene-expression
readout per tumor sample (bulk RNA-seq), you get a separate readout for every single
cell in the sample. This still matters here because averaging would hide exactly the
thing we care about: a small subpopulation of antigen-negative cells sitting inside an
otherwise antigen-positive tumor. That part of the scientific logic is untouched by
where the matrices came from.

---

## 3. The pipeline, stage by stage, in plain terms

### Stage A — Get the processed matrices (no alignment needed)
Since raw reads aren't available, there's no alignment step to run. Each sample's
Cell Ranger output — the count matrix (rows = genes, columns = cells, values = counts)
plus the barcode and feature lists that go with it — is downloaded directly from GEO's
FTP as a single ~970 MB archive, then unpacked per sample. This replaces what would
otherwise be Biowulf/SINCLAIR's job. No HPC cluster is strictly necessary for this
project anymore — the full matrix set is small enough to work with on a laptop, though
keeping it on Biowulf storage for now is fine for convenience.

One thing this stage still needs to do, that SINCLAIR would have handled automatically:
basic per-sample QC (removing low-quality cells, doublets) and integration across all
62 samples onto a common coordinate system, correcting for sample-to-sample batch
effects. This now happens locally in R (Seurat + Harmony) rather than as part of an
HPC pipeline.

### Stage B — Figure out what kind of cell each cell is
Unchanged. A bone marrow sample isn't just tumor — it's tumor cells mixed in with T
cells, NK cells, normal B cells, myeloid cells, and more. **Annotation** is the step
where we label each cell by type, based on which genes it's expressing (T cells express
CD3, plasma cells express CD38 and MZB1, etc.). We use a reference-based tool (SingleR)
plus a manual marker check as a sanity test.

### Stage C — Separate malignant plasma cells from normal ones
Unchanged. Even after identifying "plasma cells," some are the patient's cancer and
some are ordinary, healthy plasma cells still in the marrow. You can't tell them apart
by looks — you have to use **clonality**. Every antibody-producing cell commits to
making either a "kappa" or "lambda" light chain, and normal plasma cells are a healthy
mix of both. Cancer is a single runaway clone, so in an involved marrow, 90%+ of plasma
cells will share one light chain type. We use that skew to call which cells are
malignant.

### Stage D — Score each malignant cell for BCMA and GPRC5D
Unchanged. Per malignant cell, we check: does this cell express the BCMA gene? The
GPRC5D gene? Four possible categories per cell — positive for both, positive for one,
or positive for neither ("double-negative").

### Stage E — The novel metric: dual-antigen escape fraction
Unchanged. For each patient: what fraction of their malignant cells are
double-negative — i.e., would be invisible to a therapy targeting both BCMA and
GPRC5D simultaneously. High score = a dual-target CAR-T might still leave disease
behind. Low score = a dual-target approach should cover nearly the whole tumor.

### Stage F — Bring in the immune microenvironment (CellChat)
Unchanged. Split patients into high-escape-risk and low-escape-risk groups and ask,
using CellChat, whether the high-risk group also shows weaker immune cell-cell
signaling.

### Stage G — The decision packet
Unchanged. A ranked table of patients by double-negative fraction, a couple of figures,
and a short written interpretation.

---

## 4. Why this still shows real skill — revised

The original framing leaned on "starts from raw sequencing data" as the headline
differentiator. That claim no longer applies to this project, and claiming it would be
inaccurate. Here's what's still true, and worth saying instead:

- **The malignant-cell-calling step still uses biology, not just clustering** —
  light-chain restriction is a real, defensible method, not a shortcut. This is
  untouched by where the input matrices came from.
- **The core metric is still something the original dataset's authors didn't
  compute** — you're not reproducing their paper, you're asking a new question of
  their data. Also untouched.
- **It still closes the loop from biology to a business-relevant conclusion** — the
  output is directly usable for the actual strategic question (single- vs. dual- vs.
  sequential-target CAR-T design).
- **New, honest talking point**: you correctly diagnosed a real-world data
  availability problem — confirmed directly against NCBI's BioProject/SRA records
  rather than assuming — and adapted the plan accordingly instead of either giving up
  or claiming to have done something you didn't. That's a legitimate skill to name if
  asked how the project evolved, and it's a more defensible story than papering over
  the gap.
- **What's honestly no longer a differentiator**: going from FASTQ yourself. Don't
  claim this. If asked directly ("did you process raw sequencing data?"), the accurate
  answer is: "No — I checked, and this dataset's authors only deposited processed
  Cell Ranger output, not raw reads. I started from their processed matrices and did
  all the biology-specific work — QC, integration, malignant cell calling, antigen
  scoring — myself from there."

---

## 5. What to say if asked "walk me through your approach" in the interview

Revised, accurate version:
1. "I wanted to test a question adjacent to Legend's actual antigen-escape problem,
   using open myeloma scRNA-seq data."
2. "I checked whether raw sequencing data was available for reprocessing, but this
   dataset's authors only deposited processed Cell Ranger matrices — no FASTQ, no SRA
   accession. So I started from their processed per-sample matrices and did the
   biology-specific work myself downstream: QC, batch integration, malignant cell
   calling, antigen scoring."
3. "The key idea is that most antigen-escape work looks at before/after relapse. I
   instead asked how much escape risk is already present in a patient's tumor at
   baseline, before any treatment — by measuring what fraction of the malignant clone
   is negative for both BCMA and GPRC5D at once."
4. "I tied that back to the immune microenvironment with CellChat, to see whether
   high-risk patients are also immunologically colder."
5. "The output is a per-patient ranking that could inform which patients are better
   candidates for a dual-target vs. sequential-target strategy."

If pressed on why not raw FASTQ: "I looked — this study's authors didn't deposit raw
reads publicly, only processed matrices, so reprocessing from scratch wasn't an option
without contacting them directly for access."

That's the whole story, accurately stated — everything else in the technical doc is
how you'd back it up if asked for detail.
