"""Smoke tests for the ``comet`` console script.

Run in a subprocess so they exercise the real entry point, argument parsing and
exit codes rather than calling main() directly.
"""

import numpy as np
import pytest

from comet.core.io_utils import load_thunderstorm_csv
from comet.tests.conftest import run_python


def run_cli(*args):
    """Invoke the CLI via -m so it works without the console script installed."""
    return run_python(["-m", "comet.cli.cli_interface"] + list(args))


class TestArgumentHandling:
    def test_help_exits_cleanly(self):
        result = run_cli("--help")
        assert result.returncode == 0
        assert "COMET drift correction CLI" in result.stdout

    def test_missing_required_args_fails(self):
        result = run_cli()
        assert result.returncode != 0

    def test_unsupported_input_extension_fails(self, tmp_path):
        bad = tmp_path / "data.xyz"
        bad.write_text("nonsense")

        result = run_cli("-i", str(bad), "-o", str(tmp_path / "out.csv"), "-sv", "2")

        assert result.returncode != 0
        assert "Unsupported input format" in result.stderr

    def test_extension_check_is_case_insensitive(self, tmp_path, thunderstorm_csv):
        upper = tmp_path / "UPPER.CSV"
        upper.write_text(thunderstorm_csv.read_text())
        out = tmp_path / "out.csv"

        result = run_cli("-i", str(upper), "-o", str(out), "--format", "csv",
                         "-sm", "2", "-sv", "2", "-isig", "120", "-tsig", "10",
                         "-d", "100", "--mode", "cpu")

        assert result.returncode == 0, result.stderr
        assert out.exists()


class TestCsvPipeline:
    def test_end_to_end_csv_correction(self, tmp_path, thunderstorm_csv):
        out = tmp_path / "corrected.csv"

        result = run_cli(
            "-i", str(thunderstorm_csv), "-o", str(out), "--format", "csv",
            "-sm", "2", "-sv", "2", "-isig", "120", "-tsig", "10", "-d", "100",
            "--mode", "cpu",
        )

        assert result.returncode == 0, result.stderr
        assert out.exists()

        corrected = load_thunderstorm_csv(str(out))
        original = load_thunderstorm_csv(str(thunderstorm_csv))
        assert corrected.shape == original.shape
        assert np.isfinite(corrected).all()

    def test_output_differs_from_input(self, tmp_path, thunderstorm_csv):
        """A drifting dataset must actually be changed by the correction."""
        out = tmp_path / "corrected.csv"

        result = run_cli(
            "-i", str(thunderstorm_csv), "-o", str(out), "--format", "csv",
            "-sm", "2", "-sv", "2", "-isig", "120", "-tsig", "10", "-d", "100",
            "--mode", "cpu",
        )
        assert result.returncode == 0, result.stderr

        corrected = load_thunderstorm_csv(str(out))
        original = load_thunderstorm_csv(str(thunderstorm_csv))
        assert not np.allclose(corrected[:, :2], original[:, :2])


class TestCrossFormat:
    def test_h5_input_to_csv_output(self, tmp_path, drifting_dataset):
        """--format csv with an .h5 input used to hand the h5 path to read_csv."""
        from comet.core.io_utils import save_dataset_as_ms_h5

        locs, _ = drifting_dataset
        source = tmp_path / "input.h5"
        save_dataset_as_ms_h5(locs[:, :3], locs[:, 3], 100.0, filename=str(source))
        out = tmp_path / "out.csv"

        result = run_cli("-i", str(source), "-o", str(out), "--format", "csv",
                         "-sm", "2", "-sv", "2", "-isig", "120", "-tsig", "10",
                         "-d", "100", "--mode", "cpu")

        assert result.returncode == 0, result.stderr
        assert out.exists()

        corrected = load_thunderstorm_csv(str(out))
        assert len(corrected) == len(locs)
        assert np.isfinite(corrected).all()


class TestH5Pipeline:
    def test_h5_output_honours_the_output_flag(self, tmp_path, thunderstorm_csv):
        """--format h5 must write to --output, not open a file dialog."""
        out = tmp_path / "corrected.h5"

        result = run_cli(
            "-i", str(thunderstorm_csv), "-o", str(out), "--format", "h5",
            "-sm", "2", "-sv", "2", "-isig", "120", "-tsig", "10", "-d", "100",
            "--mode", "cpu",
        )

        assert result.returncode == 0, result.stderr
        assert out.exists(), "h5 output was not written to the --output path"
