"""
mm_escape — MM dual-antigen (BCMA/GPRC5D) escape risk analysis.

An importable library, NOT a pipeline sequence: modules are named by function and are
never numbered, because several (plotting, robustness, gene_space) are used by more
than one stage. The numbered things are notebooks and their `results/NN_*/` output
directories, which match 1:1. See the main project document.

Submodules are imported lazily — `gene_space` needs only pandas, while `qc` pulls the
rpy2/scDblFinder bridge and `communication` pulls LIANA+, and those live in different
conda environments. Importing the package must not require every environment at once.
"""

__version__ = "0.1.0"

__all__ = ["config", "gene_space", "io"]


def __getattr__(name: str):
    if name in __all__:
        import importlib

        return importlib.import_module(f".{name}", __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
