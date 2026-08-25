# Stage 05b — integration benchmark decision

## The rule, declared before running

A method replaces the incumbent only if **all four** hold. Thresholds are in
`benchmark.DECISION_TOLERANCES` and were fixed before any arm was run.

1. `batch_improved` — immune batch removal actually improves.
2. `bio_preserved` — broad cell identity is not sacrificed for that mixing.
3. `depth_ok` — R²(depth ~ latent) does not rise materially. R² is used
   because it depends only on the embedding's column span and is therefore
   rotation-invariant; latent axes are arbitrary across methods.
4. `overcorrection_ok` — plasma-cell mixing never earns a method anything. A
   jump alone flags it; a jump **with** rising depth association disqualifies
   it, that pair being the signature of the censoring being smoothed over.

Scoring is on the **immune compartment**. Global scIB is reported as a
secondary reference and is never used for selection.

## Verdict

**No arm qualified. `harmony_stage05` stays.** That is a real result, not a failure to find one — the incumbent survived a comparison it could have lost.

## Per-arm results

| arm | batch_score | bio_score | depth_r2 | plasma_mixing | batch_improved | bio_preserved | depth_ok | plasma_flagged | overcorrection_ok | eligible |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| unintegrated | 0.4503 | 0.6914 | 0.6603 | 0.01439 | yes | yes | no | no | yes | no |
| harmony_sample | 0.6155 | 0.7181 | 0.5093 | 0.515 | yes | yes | no | yes | no | no |
| scvi_sample | 0.5695 | 0.7013 | 0.5412 | 0.452 | yes | yes | no | yes | no | no |
| scanorama_sample | 0.4504 | 0.7228 | 0.6904 | 0.1605 | yes | yes | no | no | yes | no |
| harmony_stage05 | 0.4272 | 0.7004 | 0.3691 | 0.03823 | no | yes | yes | no | yes | no |
| harmony_cohort | 0.5909 | 0.7064 | 0.6072 | 0.7706 | yes | yes | no | yes | no | no |
| scvi_cohort | 0.4922 | 0.69 | 0.5757 | 0.01743 | yes | yes | no | no | yes | no |

## What this benchmark cannot do

**No integration method restores cells that were never deposited.** WashU
cohorts 1 and 2 were cut at 10,000 UMIs before deposit; the high-RNA portion
of their plasma-cell population is absent from the counts, whichever arm
wins. A well-mixed latent space has not undone that ascertainment bias.
**Stage 08 still owes its truncate-all-cohorts-at-10,000 sensitivity
analysis, regardless of the outcome here.**

Nor can this move the headline metric: the embedding feeds only stages 06 and
11, antigen calls read `layers['counts']`, and malignant subclustering is
per-patient and un-integrated.

Metrics computed on a stratified subsample of 40,021 cells.
