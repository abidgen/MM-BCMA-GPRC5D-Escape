# What We're Doing, and Why — Plain-Language Overview

**Revision note:** This project was originally built substantially in R (Seurat,
Harmony, CellChat), reached a working state through QC/doublet-removal on the full
sample cohort, then was restarted from scratch in Python (scanpy, harmonypy, LIANA+)
for better fit with the author's primary tooling. Every piece of hard-won knowledge
about the data itself — file-format quirks, a mixed-reference-build problem, a
patient-ID mapping gap — carries over unchanged; only the code does not. This doc's
scientific logic (Stages B through G below) is untouched by that switch.

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
~970 MB total, about 204,000 cells before quality filtering. The raw sequencing reads
exist but sit behind a controlled-access application at dbGaP, and — more importantly —
the *unfiltered* count matrices (the ones that include empty droplets) were never
deposited anywhere. That second point matters more than it might sound: a common
cleanup step called ambient RNA correction needs those unfiltered matrices, so it
simply isn't available for this dataset no matter what access you have. We handle that
gap differently instead — more on that in Stage D below.

**One wrinkle worth knowing about the cohorts.** The 62 samples were collected at
different sites on different generations of the 10x kit, and they differ by roughly
1.9× in how many genes get detected per cell. That is a problem here specifically
because the number this project reports is a *fraction of zeros* — if one cohort's
cells are read less deeply, more of their genes look absent, and that can masquerade as
biology. So which cohort a sample came from is carried through the whole analysis as a
variable to correct for, not discarded.

scRNA-seq still means what it always meant: instead of one averaged gene-expression
readout per tumor sample (bulk RNA-seq), you get a separate readout for every single
cell in the sample. This still matters here because averaging would hide exactly the
thing we care about: a small subpopulation of antigen-negative cells sitting inside an
otherwise antigen-positive tumor.

---

## 3. The pipeline, stage by stage, in plain terms

### Stage A — Get the processed matrices, clean up cell-level quality
Each sample's Cell Ranger output — the count matrix (rows = genes, columns = cells,
values = counts) plus barcode and gene lists — is downloaded from GEO and loaded.
Some cells are broken or dying (leaky membranes, mostly-mitochondrial RNA left behind)
and get filtered out; some "cells" are actually two cells stuck together in one droplet
(doublets) and get flagged and removed too. Rather than picking one fixed cutoff number
for "how broken is too broken," we use a statistical method that looks at the actual
spread of quality across all the real cells in this dataset and flags anything unusually
far from the middle of that spread — this adapts to the data instead of imposing an
arbitrary number that might not fit. That spread is measured *within each cohort*
rather than across all of them at once, since otherwise the shallower cohort's normal
cells would all look like outliers relative to the deeper cohort's.

Two things about the files themselves were not what the project believed, and both were
caught by checking rather than by assuming. The gene lists carry only gene *names*, not
stable identifiers — and gene names drift between reference versions, so matching
samples on names both loses genes and, worse, can match the *wrong* gene. The
identifiers turned out to be reconstructible from the public references the files were
built from, which recovered about 33,000 usable genes instead of 22,000. And one
sample, long excluded on the grounds that its reference was missing the BCMA gene
entirely, turned out to be a perfectly ordinary sample whose gene-name file had simply
been cut off part-way through being written; BCMA was past the cut, not absent. It is
repaired and kept.

### Stage B — Figure out what kind of cell each cell is
A bone marrow sample isn't just tumor — it's tumor cells mixed in with T cells, NK
cells, normal B cells, myeloid cells, and more. **Annotation** is the step where we
label each cell by type, based on which genes it's expressing.

There are two ways to do this, and we do both. The **manual** way is to look up which
genes are known to mark each cell type and check which clusters of cells are expressing
them. The **automatic** way is to use a tool trained on millions of already-labelled
cells and let it assign labels for you. Automatic is obviously faster — the question is
whether it's *right* on this particular dataset, and you can't know that by trusting it.

So we run the manual version, two different automatic tools, and then compare all three.
The decisive check is simple enough to eyeball: take the known marker genes and see
whether the automatic labels line up with them. If the cells an automatic tool calls
"T cells" are exactly the cells expressing the known T-cell genes, then it has already
captured everything the manual approach would have — and doing it by hand is just
labour. Where it doesn't line up, we keep the manual answer for that cell type.

Two details that matter more than they look:

- **The decision is made separately for each cell type, not once overall.** The tools
  are expected to fail on *different* cell types — the popular automatic reference is
  built from immune cells, so it handles T/NK/B/myeloid well but has never seen the
  red-blood-cell precursors and stem cells that live in marrow. Judging all-or-nothing
  would mean throwing away good labels for one population because of an unrelated
  weakness in another.
- **We write down the pass/fail bar before looking at any results.** Otherwise "pick
  the best method" quietly turns into "justify whichever one looked nicer." The bar is
  strictest for plasma cells, because those are the cells the entire project is
  measuring — get that boundary wrong and every later number inherits the error.

One thing no automatic tool can do here: tell a *cancerous* plasma cell from a healthy
one. Both look like plasma cells, because the reference data these tools learned from
only ever contained healthy ones. That's fine — it's Stage C's job, and it needs a
completely different kind of evidence.

We also record, separately from the cell-type label, what each cell is *doing* — is it
dividing, is it under stress, is it responding to interferon. A cell has one identity
but can be doing several of those at once, so they're kept as sliding scales rather than
being crushed into the identity label. A dividing plasma cell is a plasma cell that
happens to be dividing, not a third kind of cell.

### Stage C — Separate malignant plasma cells from normal ones
Unchanged. Even after identifying "plasma cells," some are the patient's cancer and
some are ordinary, healthy plasma cells still in the marrow. You can't tell them apart
by looks — you have to use **clonality**. Every antibody-producing cell commits to
making either a "kappa" or "lambda" light chain, and normal plasma cells are a healthy
mix of both. Cancer is a single runaway clone, so in an involved marrow, 90%+ of plasma
cells will share one light chain type.

### Stage D — Score each malignant cell for BCMA and GPRC5D
Mostly unchanged, with one refinement. Per malignant cell, we check: does this cell
express the BCMA gene? The GPRC5D gene? Because we can't do the usual cleanup for
background contamination (Stage A's note above — plasma cells in particular leak a lot
of antibody-gene RNA into their surroundings, which can make a truly negative cell look
faintly positive by mistake), we set the bar for "counts as positive" a little above
zero — specifically, above whatever background level shows up in cell types that have
no business expressing these genes at all (T cells, immune cells that aren't plasma
cells). That gives us a more trustworthy positive/negative call than just checking for
any signal greater than zero.

### Stage D2 — Make sure the "negative" calls are real
This is the part that got substantially strengthened in a later design review, and
it's worth understanding why. The headline number is a **fraction of zeros** — and a
zero in this kind of data can mean two completely different things. It can mean the
cell genuinely doesn't express the gene. Or it can mean the measurement simply missed
it, which happens constantly in single-cell data and is called *dropout*.

Those two failure modes push the answer in opposite directions. Background
contamination (Stage D) makes truly negative cells look faintly positive, which makes
the escape fraction look **too small**. Dropout makes truly positive cells look
negative, which makes it look **too big**. The original plan only accounted for the
first one. The second is actually the bigger worry here, because GPRC5D is a gene
that's expressed at low levels even when it *is* expressed, so it's exactly the kind
of gene the measurement misses.

Rather than pick a number and hope, the analysis now:
- computes the escape fraction under several different "counts as positive"
  thresholds, and reports whether the **patient ranking stays the same** across all of
  them (the ranking surviving is the real result, not any single number);
- checks whether patients with deeper sequencing get systematically lower escape
  fractions — if they do, the metric is secretly measuring sequencing depth rather
  than biology, and that would need to be caught *before* showing anyone the ranking;
- puts an error bar on every patient, since some patients contribute fifteen times as
  many cells as others and a bare ranking hides that completely;
- cross-checks against a **completely separate bulk RNA-seq measurement** of the same
  samples, which was downloaded at the start of the project and originally unused. If
  the bulk data shows the antigen where the single-cell data shows nothing, that's
  direct evidence of how much dropout is happening. Note what this can and can't do:
  bulk tells you how much of each antigen is around overall, but it averages every cell
  together, so it can't see *which cells* have which antigen. A tumor that's half
  BCMA-only and half GPRC5D-only looks, in bulk, like a tumor with plenty of both — even
  though not a single cell carries both. So bulk checks the antigen levels and whether
  the single-cell zeros are believable; it can't check the escape fraction itself.

### Stage E — The novel metric: dual-antigen escape fraction
For each patient: what fraction of their malignant cells are double-negative — i.e.,
would be invisible to a therapy targeting both BCMA and GPRC5D simultaneously.

### Stage E1.5 — Are the *same* cells losing both antigens?
This is a sharper question than the escape fraction itself, and it costs almost nothing
to ask.

Suppose 30% of a patient's tumor cells lack BCMA and 20% lack GPRC5D. If those two
failures are unrelated, you'd expect 0.30 × 0.20 = 6% of cells to lack both by pure
coincidence. So a 6% escape fraction in that patient is exactly what independence
predicts — two separate, partial problems that happen to overlap sometimes. Now suppose
a different patient also has 6% double-negative cells, but their individual rates only
predict 1.5%. That's four times more overlap than chance, and it means something quite
different: the *same cells* are shutting down both targets together.

What does that mean practically? Adding a second target still helps the second patient
— just less than you'd hope. If 30% of cells lack BCMA, adding GPRC5D brings the
uncovered share down to whatever fraction lacks *both*: 6% if the two failures are
unrelated, maybe 15% if they travel together. Fifteen percent is worse than six, but
it's still half of thirty. So co-escape doesn't measure whether the second target is
worth adding — it measures **how much of the benefit you expected gets eaten by the two
failures overlapping**. (An earlier draft of this document said the second patient is
"the one dual targeting doesn't help." That was too strong, and the arithmetic above is
why.)

So each patient gets two numbers here, kept separate because they answer different
questions: a **co-escape enrichment** figure (observed double-negatives divided by what
chance alone predicts — a fact about biology) and an **incremental coverage gain**
(how much the second target actually reduces the uncovered fraction — a fact about
clinical value). A patient can score high on both.

One trap here, which is worth knowing about because it points the wrong way. Cells that
were sequenced shallowly show false zeros for *both* genes at once, purely because
there wasn't enough data. That alone creates fake "co-escape" — the same cells reading
negative for both, for a purely technical reason. So the comparison is done within
groups of similarly-sequenced cells, and the uncorrected number is reported next to the
corrected one so the size of that artifact is visible rather than hidden.

### Stage E2 — Is the escape population *structured*, or just noise?
**This is arguably the real scientific question.**

Saying "3% of this patient's tumor cells are double-negative" and saying "this patient
already has a 3% resistant subclone" sound identical but aren't. Only the second one
means anything clinically — because a *subclone* is a coherent group of related cells
that therapy would actively select for, leaving them to regrow. Cells that just
happen to be scattered randomly across the tumor are far more likely to be measurement
noise, and don't predict relapse the same way.

The tempting shortcut is: if the double-negative cells cluster together in "cell
similarity space", it's a subclone; if they're sprinkled at random, it's dropout. That
shortcut is too quick, and the plan used to rely on it. Cells can resemble each other
for lots of reasons that have nothing to do with being genetically related — they might
all be dividing, all stressed, all responding to interferon, or all just sequenced
shallowly. And a genuinely related group of cells doesn't have to look alike.

So the claim is built in three steps, each one stronger than the last:
1. **Are the double-negative cells non-randomly placed?** If yes, something structures
   them — that alone rules out pure scatter.
2. **Do they share a distinctive gene program?** If yes, there's an escape-associated
   *state* — a recognizable mode these cells are in.
3. **Do they share genetic changes** (chromosome gains and losses, read off the
   expression data)? Only if this holds does the word *subclone* get used.

The honest caveat: step 3 is hard. Spotting the difference between a tumor cell and a
healthy cell is much easier than spotting differences *within* one patient's tumor, and
at this data's depth it often won't be resolvable at all. When that happens the answer
is "couldn't tell", not "no subclone" — because quietly treating an underpowered test
as a negative result would systematically understate exactly the risk this project
exists to measure.

We also look at *what else* is different about those escape cells — what genes they
turn up or down. One specific hypothesis is written down in advance (so it stays a
real test and not a story invented after seeing the answer): an enzyme called
γ-secretase physically cuts BCMA off the surface of myeloma cells, and drugs that
block it are already being trialed alongside BCMA CAR-T for exactly that reason. If
the escape cells turn out to be high in the γ-secretase machinery, that points at a
mechanism *and* an existing intervention.

### Stage F — Bring in the immune microenvironment
Changed in tool and in statistics. We use LIANA+, a Python tool that includes the same
underlying method as CellChat (the R tool originally planned) plus several others, so
the comparison is if anything more robust than relying on one method alone. The
question: do high-escape-risk patients also show weaker immune cell-cell signaling?

The original plan was to split patients into a high group and a low group and compare
the two piles of cells. That turns out to be a statistical trap called
**pseudoreplication** — pooling cells across patients treats one patient's few thousand
cells as a few thousand independent data points, when really they're all just *one
patient*. It makes almost anything look significant. So instead each patient gets one
score, the comparison is across ~41 patients rather than ~200,000 cells, and escape
risk is used as a sliding scale instead of being chopped into high/low.

There's also an obvious alternative explanation to rule out: maybe high-escape patients
simply *have* fewer T and NK cells to begin with, in which case "weaker immune
signaling" is just an artifact of there being fewer immune cells to signal. So immune
cell abundance gets measured separately and controlled for.

### Stage F2 — What normal marrow cells carry these antigens
The healthy bone marrow samples in the dataset were originally going to sit unused.
They're now doing two jobs.

First, they're a **sanity check on the whole method**. Healthy marrow has no cancer
clone in it, so when we run the malignant-cell-finding logic on healthy samples, it
should find nothing. If it "discovers" a tumor in a healthy person, the method is
broken and every other number is worthless. Cheap test, and it either passes or it
saves the project.

Second, they tell us whether *normal* plasma cells also carry BCMA and GPRC5D — which
matters because a CAR-T can't tell the difference. A target that's on the tumor but
also on healthy tissue causes side effects. This is a real clinical distinction between
the two antigens: BCMA sits on normal plasma cells and B cells fairly broadly, while
GPRC5D is more tumor-specific in marrow but shows up in skin and nails, which is where
the known nail and rash side effects of GPRC5D-targeted drugs come from. Adding this
turns the project from "which target kills more tumor" into "which target kills more
tumor *per unit of harm*," which is the actual question.

### Stage F3 — Would a different pair of targets do better?
Rather than only asking about BCMA and GPRC5D, we score several other candidate myeloma
targets too and compute, for every possible pair or trio, how much of each patient's
tumor would be left uncovered. That answers something the two-antigen number
structurally can't: *for this specific patient, is BCMA + GPRC5D even the right
combination, or would a different pairing cover more of their disease?* And coverage
gets weighed against the normal-tissue expression from Stage F2 — covering 100% of the
tumor isn't a win if the target is also all over healthy tissue.

### Stage G — The decision packet
A ranked table of patients by double-negative fraction — now with **error bars**, since
a plain bar chart would claim a precision this measurement doesn't have — plus the
target-coverage comparison, the "is it a subclone" score, the bulk-data cross-check, a
few figures, and a short written interpretation.

It also states plainly what the analysis *can't* say. The biggest one: CAR-T therapy
binds to **protein on the cell surface**, and this analysis measures **RNA inside the
cell**. Those usually track each other, but not perfectly — and for BCMA specifically
there's a known mechanism that actively strips the protein off the surface while the
RNA stays put. That's the first thing a knowledgeable person will ask about, so it goes
in the deliverable rather than waiting to be caught.

---

## 4. Why this still shows real skill

- **The malignant-cell-calling step still uses biology, not just clustering** —
  light-chain restriction is a real, defensible method, not a shortcut.
- **The core metric is still something the original dataset's authors didn't
  compute** — this asks a new question of their data, not a reproduction of their paper.
- **It still closes the loop from biology to a business-relevant conclusion.**
- **Recognizing and working around real data limitations** — raw reads only under
  dbGaP controlled access, no unfiltered matrices for ambient correction anywhere, a
  mixed-reference-build problem, three collection cohorts on different 10x chemistries,
  and a patient-ID mapping gap — documenting each one honestly with a stated
  mitigation, rather than either ignoring the gap or claiming a workaround fixes it
  completely.
- **Reading the deposit sceptically, and being repaid for it.** Two of this project's
  load-bearing "facts" were wrong. The gene files carry no Ensembl IDs, but they turned
  out to be positional dumps of public references, so the IDs were reconstructible —
  recovering 32,991 genes against 22,164 on symbols, and preventing a silent
  mis-pairing where the same symbol means different genes in the two builds. And one
  sample, excluded for years as an incompatible reference missing BCMA, was really a
  normal sample whose gene file had been truncated mid-write; the "missing" genes were
  past the cut. Both were found by checking a claim against the files instead of
  inheriting it.
- **Auditing your own headline metric and finding it under-defended.** The original
  plan accounted for one source of error and missed the larger, opposite-signed one
  (dropout). Catching that in your own design — and then reporting a range and a
  ranking-stability claim rather than a single confident number — is a more useful
  signal about how someone works than any individual result.
- **Distinguishing "3% of cells are negative" from "there's a 3% resistant subclone."**
  These are easy to conflate and only one of them predicts relapse. Knowing they're
  different — and then *not* over-claiming the second when the data only supports the
  first — is the difference between a descriptive number and an actual finding.
- **Asking whether the same cells lose both antigens**, rather than stopping at how
  many lose both. It's one contingency table per patient, and it separates two tumors
  that look identical on the headline number. Note what it does *not* say: correlated
  loss erodes how much complementarity the second target delivers, it does not make the
  second target worthless — adding GPRC5D to BCMA still moves the uncovered fraction
  from `P(BCMA⁻)` to `P(BCMA⁻ ∩ GPRC5D⁻)`, which is a real gain even under strong
  co-loss. The incremental coverage gain is reported next to the enrichment for exactly
  that reason.
- **Refusing to correct what can only be bounded.** The obvious dropout "correction"
  multiplies each cell's two negativity probabilities — which assumes the independence
  the co-escape test exists to interrogate. It is reported as a technical baseline to
  compare against, never as a corrected truth, and no dropout-corrected point estimate
  is claimed.
- **Using the unglamorous data.** The matched bulk RNA-seq and the healthy-marrow
  controls were both already downloaded and both originally slated to sit unused. They
  turned out to supply the project's only independent check on its antigen levels and
  its only measurement of what normal marrow plasma cells express.
- **Getting the statistics right on the immune comparison** — treating patients rather
  than cells as the unit of replication, and controlling for the obvious alternative
  explanation instead of hoping nobody raises it.
- **Rebuilding in the tool you're actually fluent in, mid-project, rather than pushing
  through in a less comfortable one** is itself a reasonable engineering call worth being
  able to explain plainly if asked — comfort with your tools measurably affects how many
  bugs you catch versus miss, and that's a real factor, not an excuse.

---

## 5. What to say if asked "walk me through your approach" in an interview

1. "I wanted to test a question adjacent to Legend's actual antigen-escape problem,
   using open myeloma scRNA-seq data."
2. "I checked whether raw sequencing data was available, but this dataset's authors
   only deposited processed, filtered matrices — no FASTQ, no unfiltered counts either.
   That second gap meant I couldn't run standard ambient-RNA cleanup, so I built a
   background-noise-derived detection threshold instead, rather than ignoring the issue."
3. "The key idea is that most antigen-escape work looks at before/after relapse. I
   instead asked how much escape risk is already present in a patient's tumor at
   baseline, before any treatment — by measuring what fraction of the malignant clone
   is negative for both BCMA and GPRC5D at once."
4. "The number I'm computing is a fraction of zeros, which is the most fragile thing
   you can measure in single-cell data — a zero is either real biology or a missed
   measurement. So I bounded both error directions instead of one, reported the
   ranking as an interval, and validated the antigen calls against matched bulk
   RNA-seq from the same samples — being careful to say that bulk checks the antigen
   *levels*, not the escape fraction, since averaging every cell together destroys the
   per-cell combination that the whole metric is about."
5. "Then the two parts I actually care about. First: are the *same* cells losing both
   antigens more often than chance predicts from their individual rates? Because if the
   two antigens fail independently, dual targeting works — and if they fail together in
   the same cells, it doesn't. Second: is that escape population structured or scattered,
   which I tested in three escalating steps and only called a subclone at the step that
   involves shared genetic changes. Transcriptional similarity alone has too many
   innocent explanations to carry that word."
6. "I also used the healthy marrow controls to check whether normal plasma cells carry
   these antigens too. That's marrow expression context, not a safety profile — GPRC5D's
   real clinical toxicity is in skin and nails, and a bone marrow dataset can't see that
   at all."
7. "I tied that back to the immune microenvironment as an explicitly exploratory
   extension — whether high-risk patients are also immunologically colder, with the
   patient as the unit of replication and immune-cell abundance controlled. Forty-one
   patients against hundreds of possible signaling pairs doesn't support a confirmatory
   claim, so I don't make one."
8. "The output is per-patient risk *tiers* rather than a ranked list — the error bars
   overlap too much for '#1, #2, #3' to be honest — with the co-escape and structure
   results as separate columns, plus a coverage comparison across other candidate
   targets in case a different pairing covers a given patient's disease better."

---

## 6. Planned next phase — checking whether the finding replicates

Once the primary analysis is fully finished, the plan is to repeat the same logic on a
second, independent dataset (a different research group, a different lab technology) to
see whether the core finding shows up there too. This is deliberately sequenced *after*
finishing the primary analysis — a complete single-dataset analysis is worth more than a
rushed two-dataset one.
