"""Guards on the installed package itself.

These are the checks that catch a broken *distribution* rather than broken
maths: an unimportable module, a missing export, an undeclared dependency, or a
version that disagrees with the installed metadata.
"""

import importlib
import pkgutil
import sys

import pytest

import comet
from comet.tests.conftest import run_python

# Every module must import on a machine with no GPU and no optional extras.
PUBLIC_MODULES = [
    "comet",
    "comet.core.backends",
    "comet.core.cpu_wrapper",
    "comet.core.drift_optimizer",
    "comet.core.experimental",
    "comet.core.interpolation",
    "comet.core.io_utils",
    "comet.core.pair_indices",
    "comet.core.qc_utils",
    "comet.core.segmenter",
    "comet.core._dialogs",
    "comet.cli.cli_interface",
    "comet.utilities.post_analysis_utilities",
    "comet.batch.mfx_batch",
    "comet.tests.self_test",
]


class TestImports:
    @pytest.mark.parametrize("module", PUBLIC_MODULES)
    def test_module_imports(self, module):
        importlib.import_module(module)

    def test_every_submodule_is_importable(self):
        """Walk the package so a newly added module cannot silently break."""
        failures = []
        for info in pkgutil.walk_packages(comet.__path__, prefix="comet."):
            # cuda/torch wrappers legitimately need their optional backend
            if "cuda_wrapper" in info.name or "pytorch_wrapper" in info.name:
                continue
            try:
                importlib.import_module(info.name)
            except Exception as exc:
                failures.append("{}: {}".format(info.name, exc))
        assert not failures, "modules failed to import: {}".format(failures)

    def test_import_does_not_require_tkinter(self):
        """tkinter is missing from many Linux and container Pythons."""
        code = (
            "import sys\n"
            "sys.modules['tkinter'] = None\n"   # force any tkinter import to fail
            "import comet\n"
            "from comet.core import io_utils, drift_optimizer\n"
            "print('ok')\n"
        )
        result = run_python(["-c", code], timeout=300)
        assert result.returncode == 0, result.stderr
        assert "ok" in result.stdout

    def test_dialog_helper_reports_missing_tkinter_clearly(self):
        code = (
            "import sys\n"
            "sys.modules['tkinter'] = None\n"
            "from comet.core._dialogs import ask_open_filename\n"
            "try:\n"
            "    ask_open_filename()\n"
            "except RuntimeError as e:\n"
            "    assert 'tkinter' in str(e), e\n"
            "    print('ok')\n"
        )
        result = run_python(["-c", code], timeout=300)
        assert result.returncode == 0, result.stderr
        assert "ok" in result.stdout

    def test_no_undeclared_third_party_imports(self):
        """Catch stray imports like the `from sympy import false` that once shipped.

        Scans the source statically: inspecting live module namespaces picks up
        whatever pytest's assertion rewriting injects.
        """
        import ast
        import importlib.util
        import pathlib

        declared = {
            "numpy", "scipy", "numba", "pandas", "matplotlib", "h5py",  # runtime deps
            "torch",                                                    # optional extra
            "pytest",                                                   # test extra
            "comet",                                                    # ourselves
        }

        package_root = pathlib.Path(comet.__file__).parent
        offenders = []
        for source in sorted(package_root.rglob("*.py")):
            # Python source is UTF-8 by definition (PEP 3120). Reading with the
            # platform default instead broke on Windows: cp1252 cannot decode
            # the Greek letters and arrows in the original docstrings.
            # utf-8-sig additionally tolerates a BOM, which plain utf-8 would
            # hand to ast.parse as a stray character.
            tree = ast.parse(source.read_text(encoding="utf-8-sig"), filename=str(source))
            imported = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported.update(alias.name.split(".")[0] for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                    imported.add(node.module.split(".")[0])

            for name in sorted(imported - declared):
                spec = importlib.util.find_spec(name) if name not in sys.builtin_module_names else None
                origin = getattr(spec, "origin", None) or ""
                # third-party packages live in site-packages; stdlib does not
                if "site-packages" in origin or "dist-packages" in origin:
                    offenders.append("{} imports undeclared '{}'".format(
                        source.relative_to(package_root), name))

        assert not offenders, offenders


class TestPublicApi:
    def test_version_is_a_string(self):
        assert isinstance(comet.__version__, str)
        assert comet.__version__.count(".") >= 1

    def test_version_matches_installed_metadata(self):
        try:
            from importlib.metadata import PackageNotFoundError, version
        except ImportError:
            pytest.skip("importlib.metadata needs Python 3.8+")
        try:
            installed = version("py-comet")
        except PackageNotFoundError:
            pytest.skip("py-comet is not installed in this environment")
        assert installed == comet.__version__, (
            "installed metadata says {!r} but comet.__version__ is {!r}. If you have "
            "switched release lanes or run an editable install, a stale "
            "Python_interface/py_comet.egg-info may be shadowing the real metadata; "
            "delete it and reinstall.".format(installed, comet.__version__)
        )

    @pytest.mark.parametrize("name", [
        "comet_run_kd",
        "segmentation_wrapper",
        "load_thunderstorm_csv",
        "load_normal_molecule_set",
        "save_dataset_as_thunderstorm_csv",
        "save_dataset_as_ms_h5",
        "correct_and_save_thunderstorm_csv",
        "available_backends",
        "best_backend",
        "cuda_available",
        "torch_available",
        "describe_backends",
    ])
    def test_export_is_present_and_callable(self, name):
        assert name in comet.__all__
        assert callable(getattr(comet, name))

    def test_all_entries_actually_exist(self):
        missing = [name for name in comet.__all__ if not hasattr(comet, name)]
        assert not missing, "__all__ lists non-existent attributes: {}".format(missing)


class TestBackendDetection:
    def test_cpu_is_always_available(self):
        assert "cpu" in comet.available_backends()

    def test_best_backend_is_one_of_the_available(self):
        assert comet.best_backend() in comet.available_backends()

    def test_probes_return_booleans_and_never_raise(self):
        assert isinstance(comet.cuda_available(), bool)
        assert isinstance(comet.torch_available(), bool)

    def test_describe_backends_mentions_all_three(self):
        text = "\n".join(comet.describe_backends())
        for backend in ("cuda", "torch", "cpu"):
            assert backend in text

    def test_mps_torch_is_never_auto_selected(self):
        """MPS is a real GPU but is slower than the compiled CPU kernel here.

        Measured across realistic segment counts: 3.4 vs 23.8 ms at 496 k pairs,
        19.7 vs 108.6 ms at 4 M, 57.0 vs 349.5 ms at 12.5 M. Eager torch
        materialises every intermediate over all pairs, so it is memory-bound.

        The invariant is that torch without CUDA is never the automatic choice
        -- not that the choice is "cpu", which would wrongly fail on a machine
        that also has a CUDA device.
        """
        from comet.core.backends import torch_device

        if not comet.torch_available() or torch_device() != "mps":
            pytest.skip("this machine has no MPS torch device")
        assert comet.best_backend() != "torch"

    def test_cpu_only_torch_is_never_auto_selected(self):
        """A CPU-device torch install must not be picked over the numba kernel.

        Asserting best_backend() == "cpu" here was wrong: on a machine with
        both CUDA and a CPU-only torch (a CUDA workstation without a CUDA
        torch build), the correct answer is "cuda". The invariant is only
        that it is never "torch".
        """
        from comet.core.backends import torch_device

        if not comet.torch_available() or torch_device() != "cpu":
            pytest.skip("torch here is not CPU-only")
        assert comet.best_backend() != "torch"


class TestConsoleScripts:
    def test_self_test_passes(self):
        """The post-install smoke test must succeed on whatever backend exists."""
        result = run_python(["-m", "comet.tests.self_test", "--mode", "cpu"])
        assert result.returncode == 0, result.stdout + result.stderr
        assert "PASSED" in result.stdout

    def test_self_test_reports_failure_with_nonzero_exit(self):
        """Forcing an unavailable backend must fail loudly, not silently pass."""
        if comet.cuda_available():
            pytest.skip("this machine has CUDA, so the cuda backend will succeed")
        result = run_python(["-m", "comet.tests.self_test", "--mode", "cuda"])
        assert result.returncode != 0
        assert "FAILED" in result.stdout
