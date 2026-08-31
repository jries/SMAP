# Contributing to COMET

Thanks for your interest in COMET. This document covers the Python package in
`Python_interface/`; the CUDA/C++ sources and the IDL interface live at the
repository root and are built with CMake.

## Development setup

```bash
git clone https://github.com/gpufit/Comet
cd Comet/Python_interface
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

Check that the install works:

```bash
comet_self_test
```

This picks the fastest backend available on your machine and needs no GPU.

## Running the tests

```bash
pytest
```

The suite runs on the CPU backend and takes about 80 seconds. No GPU is
required: tests that need one are marked and skip automatically.

| Command                        | What it runs                                  |
| ------------------------------ | --------------------------------------------- |
| `pytest`                       | Everything available on this machine           |
| `pytest -m cuda`               | CUDA backend only (skipped without an NVIDIA GPU) |
| `pytest -m torch`              | PyTorch backend only (needs the `torch` extra) |
| `pytest --run-slow`            | Also the slow regression tests                 |
| `pytest --pyargs comet.tests`  | The suite as installed, from any directory     |

If you have a GPU, please run `pytest -m cuda` before opening a pull request.
CI has no GPU, so the CUDA backend is only ever exercised by contributors.

## Release lanes

COMET publishes from two branches:

| Branch          | Version | `requires-python` | Purpose                              |
| --------------- | ------- | ----------------- | ------------------------------------ |
| `master`        | 1.1.x   | `>=3.9`           | Active development                   |
| `maint/1.0.x`   | 1.0.x   | `>=3.6`           | Keeps pip from blocking 3.6/3.7      |

Both lanes share the same source; only the packaging metadata differs.
`pip install py-comet` resolves to the right one automatically, because pip
skips releases whose `requires-python` excludes the running interpreter.

The `1.0.x` lane exists so that users on an old interpreter are not refused at
install time. 3.6/3.7 are **not a supported configuration** — the dependency
stack is difficult to get working there and CI only covers it best-effort.

This means **the package source must stay Python 3.6-compatible**, so that
fixes can be backported without rewriting them. CI enforces this:

```bash
vermin --backport dataclasses --eval-annotations -t=3.6 comet/
```

No walrus operators, no f-string `=` specifiers, no `dict | dict`, no
positional-only parameters. Test files are exempt.

Only deliberate backports land on `maint/1.0.x`; new features go to `master`.

## Adding a dependency

Runtime dependencies go in `[project.dependencies]` in `pyproject.toml`, with a
lower bound only — please don't add upper caps unless there is a known
incompatibility, since they cause resolver conflicts in users' environments.

`test_package.py::test_no_undeclared_third_party_imports` statically scans the
package for imports that are not declared, and will fail if you forget.

Anything optional belongs in an extra, and must be imported lazily inside the
function that needs it so `import comet` keeps working without it.

## Code conventions

- The library must be usable non-interactively: no `input()` prompts, no
  unconditional `plt.show()`, and no `print()` outside a `display`/`verbose`
  guard. Use the `_log()` helper in `drift_optimizer.py`.
- `tkinter` is not available in every Python installation. Import file dialogs
  through `comet.core._dialogs`, which loads it on demand.
- New public functions belong in `comet/__init__.py`'s `__all__`.

## Pull requests

1. Branch from `master`.
2. Add tests covering the change.
3. Make sure `pytest` passes and `comet_self_test` still succeeds.
4. Add an entry to `CHANGELOG.md` under `## [Unreleased]`.
