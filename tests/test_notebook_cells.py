"""Every jupytext notebook cell must compile on its own.

WHY THIS EXISTS
---------------
A `.py` notebook source is ONE compilation unit; the `.ipynb` it produces is MANY.
`jupytext` treats any `# %%` as a cell boundary — including an indented one — so this:

    # %%
    if condition:
        ...
        # %%          <- indented marker, still a cell boundary
        ...

splits an `if` header away from its body. The script compiles fine, `python -m
py_compile` passes, and the notebook fails at run time with `IndentationError`.

That is not hypothetical: `notebooks/06_annotation.py` carried exactly this defect from
the moment a `if not REUSE:` reuse branch was introduced. Every stage-06 iteration ran
green because the `.py` was being executed directly as a script, while the committed
`.ipynb` — the artifact the project treats as reviewable — had been broken the whole
time. Script-only execution cannot detect it.

So the invariant is per CELL, not per file:

    every code cell of every notebooks/*.py must compile independently

This is a fast syntax/structure check. It never executes a cell.
"""

from __future__ import annotations

from pathlib import Path

import pytest

NOTEBOOK_DIR = Path(__file__).resolve().parents[1] / "notebooks"


def _notebook_sources() -> list[Path]:
    return sorted(p for p in NOTEBOOK_DIR.glob("*.py") if "# %%" in p.read_text())


def _code_cells(path: Path) -> list[tuple[int, str]]:
    """Return (index, source) for code cells, using jupytext's own parser.

    jupytext is what actually produces the `.ipynb`, so parsing with anything else
    would be testing a different thing from the one that breaks.
    """
    import jupytext

    nb = jupytext.read(path, fmt="py:percent")
    return [
        (i, cell["source"])
        for i, cell in enumerate(nb.cells)
        if cell["cell_type"] == "code" and cell["source"].strip()
    ]


def _compile_failure(path: Path, index: int, source: str) -> str | None:
    """Compile one cell; return a formatted message on failure, else None."""
    try:
        compile(source, f"{path.name}[cell {index}]", "exec")
    except SyntaxError as exc:            # IndentationError is a subclass
        lineno = exc.lineno or 1
        lines = source.splitlines()
        lo, hi = max(0, lineno - 3), min(len(lines), lineno + 2)
        context = "\n".join(
            f"    {'>>' if n == lineno else '  '} {n:>3}| {lines[n - 1]}"
            for n in range(lo + 1, hi + 1)
        )
        return (
            f"\n{path.name}: code cell {index} does not compile on its own.\n"
            f"  {type(exc).__name__}: {exc.msg} (cell line {lineno})\n"
            f"  cell source around the failure:\n{context}\n"
            f"  The script as a whole may compile; a notebook cell is its own\n"
            f"  compilation unit. Check for a `# %%` marker inside an indented block."
        )
    return None


@pytest.mark.parametrize("path", _notebook_sources(), ids=lambda p: p.name)
def test_every_notebook_cell_compiles_independently(path: Path) -> None:
    failures = [
        msg
        for index, source in _code_cells(path)
        if (msg := _compile_failure(path, index, source)) is not None
    ]
    assert not failures, "".join(failures)


def test_no_cell_marker_sits_inside_an_indented_block() -> None:
    """The specific shape that caused the defect, checked directly and cheaply.

    A `# %%` marker is only ever valid at column 0. This catches the problem even in a
    file where the resulting cells happen to compile anyway.
    """
    offenders = []
    for path in _notebook_sources():
        for n, line in enumerate(path.read_text().splitlines(), start=1):
            if line.lstrip().startswith("# %%") and line != line.lstrip():
                offenders.append(f"{path.name}:{n}: indented cell marker -> {line!r}")
    assert not offenders, "\n".join(offenders)


def test_notebook_sources_were_actually_found() -> None:
    """A silent zero-notebook glob would make every test above vacuously pass."""
    found = _notebook_sources()
    assert len(found) >= 7, f"expected the stage 01-06 notebooks, found {len(found)}"
    assert any(p.name == "06_annotation.py" for p in found)


# --------------------------------------------------------------------------
# engineered fixtures — the defect, and a valid control
# --------------------------------------------------------------------------

BROKEN = '''# %%
if True:
    x = 1
    # %%
    y = 2
'''

VALID = '''# %%
if True:
    x = 1
    y = 2

# %%
z = 3
'''


def _cells_from_text(tmp_path: Path, text: str, name: str) -> list[tuple[int, str]]:
    p = tmp_path / name
    p.write_text(text)
    return _code_cells(p)


def test_engineered_broken_notebook_is_rejected(tmp_path: Path) -> None:
    """The indented-marker shape must fail even though the combined source is readable
    Python — this is the regression case."""
    cells = _cells_from_text(tmp_path, BROKEN, "broken_nb.py")
    assert len(cells) >= 2, "jupytext should split at the indented marker"
    failures = [
        _compile_failure(tmp_path / "broken_nb.py", i, s)
        for i, s in cells
        if _compile_failure(tmp_path / "broken_nb.py", i, s)
    ]
    assert failures, "the indented-marker notebook should NOT compile cell-by-cell"
    assert "IndentationError" in failures[0]


def test_engineered_valid_notebook_passes(tmp_path: Path) -> None:
    cells = _cells_from_text(tmp_path, VALID, "valid_nb.py")
    assert len(cells) == 2
    for i, s in cells:
        assert _compile_failure(tmp_path / "valid_nb.py", i, s) is None


def test_combined_source_compiles_even_when_cells_do_not(tmp_path: Path) -> None:
    """The point of the whole file: `py_compile` on the script would NOT catch this."""
    compile(BROKEN, "broken_nb.py", "exec")          # script-level: fine
    cells = _cells_from_text(tmp_path, BROKEN, "broken2.py")
    assert any(_compile_failure(tmp_path / "broken2.py", i, s) for i, s in cells)
