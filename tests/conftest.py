"""
Shared fixtures.

TWO TIERS, and the split is deliberate. Most of this suite runs on a fresh clone with
no data at all, because the things most worth protecting — the Ensembl-ID join, the
`make.unique` reimplementation, the truncation repair, the required-gene assertions —
are exercised entirely by what is committed under `resources/`. Only the tests that
genuinely need count matrices are gated on `raw/`, and they skip rather than fail so
the suite stays runnable anywhere.

Run everything:            pytest
Skip the slow cohort pass: pytest -m "not slow"
"""

from __future__ import annotations

import sys
from pathlib import Path

import importlib.util

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from mm_escape import config  # noqa: E402

#: The four samples that cover this deposit's failure modes. Every data-backed test
#: runs over these rather than an arbitrary pick:
#:   MMRF_1695  33538 build, MMRF cohort, 10x v3.3
#:   27522_1    33694 build, WashU 1, 10x v2 (legacy symbols: WHSC1 not NSD2)
#:   BM4        33538 build, normal donor
#:   56203_1    truncated genes.tsv, repaired on read
CANONICAL_SAMPLES = ("MMRF_1695", "27522_1", "BM4", "56203_1")

RAW_AVAILABLE = config.MANIFEST_CSV.exists() and config.SAMPLES_DIR.is_dir()

requires_data = pytest.mark.skipif(
    not RAW_AVAILABLE,
    reason=f"needs the extracted deposit at {config.SAMPLES_DIR} (see notebooks/01)",
)

#: R lives only in `mm-qc` (scDblFinder) and `mm-annotation` (SingleR), on purpose —
#: see the project's env split. The suite's home is `mm-core`, which carries no R at
#: all, so R-backed tests skip there rather than failing.
#:
#: To actually exercise them you need pytest inside an R env. `pytest` is declared in
#: `envs/env-qc.yml` but the currently-built `mm-qc` predates that line, so either
#: rebuild it or rely on `notebooks/04_qc.ipynb`, which runs the bridge over all 62
#: samples and is the stronger check anyway.
requires_r = pytest.mark.skipif(
    importlib.util.find_spec("anndata2ri") is None,
    reason="needs the mm-qc environment (rpy2 + anndata2ri + scDblFinder)",
)

#: Supplementary Table S1 is a separate gate: it is a small xlsx from the journal, not
#: part of the GEO deposit, so `raw/` can be present without it. The *parsed* tables
#: are committed under `resources/`, so only tests that re-run the parser need this.
requires_s1 = pytest.mark.skipif(
    not config.S1_XLSX.exists(),
    reason=f"needs Supplementary Table S1 at {config.S1_XLSX}",
)


@pytest.fixture(scope="session")
def manifest():
    if not RAW_AVAILABLE:
        pytest.skip("no extracted deposit")
    from mm_escape import io

    return io.load_manifest()


@pytest.fixture(scope="session")
def samples(manifest):
    """The four canonical samples, loaded once for the whole session."""
    from mm_escape import io

    return {name: io.read_sample(name, manifest=manifest) for name in CANONICAL_SAMPLES}
