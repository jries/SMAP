"""Round-trip tests for the localization IO helpers.

Fixtures are generated on the fly rather than read from a checked-in data
directory, so the suite has no external data dependency.
"""

import numpy as np
import pytest

from comet.core.io_utils import (
    load_normal_molecule_set,
    load_thunderstorm_csv,
    save_dataset_as_ms_h5,
    save_dataset_as_thunderstorm_csv,
    correct_and_save_thunderstorm_csv,
)


class TestThunderstormCsv:
    def test_load_returns_xyz_frame_columns(self, thunderstorm_csv):
        locs = load_thunderstorm_csv(str(thunderstorm_csv))

        assert isinstance(locs, np.ndarray)
        assert locs.shape[1] == 4  # x, y, z, frame
        assert np.isfinite(locs).all()

    def test_load_preserves_frame_numbers(self, thunderstorm_csv, drifting_dataset):
        expected, _ = drifting_dataset
        locs = load_thunderstorm_csv(str(thunderstorm_csv))

        np.testing.assert_array_equal(locs[:, 3].astype(int), expected[:, 3].astype(int))

    def test_load_preserves_coordinates(self, thunderstorm_csv, drifting_dataset):
        expected, _ = drifting_dataset
        locs = load_thunderstorm_csv(str(thunderstorm_csv))

        # the fixture writes 3 decimal places
        np.testing.assert_allclose(locs[:, :3], expected[:, :3], atol=1e-3)

    def test_return_essentials_false_gives_dataframe(self, thunderstorm_csv):
        df = load_thunderstorm_csv(str(thunderstorm_csv), return_essentials=False)

        assert hasattr(df, "columns")
        assert "frame" in df.columns

    def test_save_then_load_round_trip(self, tmp_path, drifting_dataset):
        locs, _ = drifting_dataset
        out = tmp_path / "round_trip.csv"

        save_dataset_as_thunderstorm_csv(locs, str(out))
        assert out.exists()

        reloaded = load_thunderstorm_csv(str(out))
        np.testing.assert_allclose(reloaded[:, :3], locs[:, :3], atol=1e-3)
        np.testing.assert_array_equal(reloaded[:, 3].astype(int), locs[:, 3].astype(int))

    def test_2d_dataset_omits_z_column(self, tmp_path, drifting_dataset):
        locs, _ = drifting_dataset
        locs = locs.copy()
        locs[:, 2] = 0.0
        out = tmp_path / "two_d.csv"

        save_dataset_as_thunderstorm_csv(locs, str(out))

        header = out.read_text().splitlines()[0]
        assert "z [nm]" not in header

        reloaded = load_thunderstorm_csv(str(out))
        assert np.all(reloaded[:, 2] == 0.0)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_thunderstorm_csv(str(tmp_path / "does_not_exist.csv"))


class TestMoleculeSetH5:
    def test_save_then_load_round_trip(self, tmp_path, drifting_dataset):
        locs, _ = drifting_dataset
        out = tmp_path / "round_trip.h5"
        pixelsize_nm = 100.0

        save_dataset_as_ms_h5(locs[:, :3], locs[:, 3], pixelsize_nm, filename=str(out))
        assert out.exists()

        reloaded = load_normal_molecule_set(str(out))
        assert reloaded.shape == (len(locs), 4)
        np.testing.assert_array_equal(reloaded[:, 3].astype(int), locs[:, 3].astype(int))

    def test_coordinates_survive_the_round_trip(self, tmp_path):
        """x must stay x. The writer used to swap the first two columns."""
        out = tmp_path / "axes.h5"
        # deliberately distinct per axis so a swap cannot go unnoticed
        coords = np.array([[100.0, 500.0, 900.0],
                           [200.0, 600.0, 800.0],
                           [300.0, 700.0, 700.0]])
        frames = np.array([0, 1, 2])

        save_dataset_as_ms_h5(coords, frames, 100.0, filename=str(out))
        reloaded = load_normal_molecule_set(str(out))

        # float32 storage in pixel units, so allow a small tolerance
        np.testing.assert_allclose(reloaded[:, :3], coords, rtol=1e-5)

    def test_separate_z_pixel_size_survives_the_round_trip(self, tmp_path):
        """z is stored in z-pixel units, so it must be scaled by the z pixel size."""
        out = tmp_path / "aniso.h5"
        coords = np.array([[100.0, 500.0, 900.0],
                           [200.0, 600.0, 800.0]])
        frames = np.array([0, 1])

        save_dataset_as_ms_h5(coords, frames, 100.0, pixelsize_z_nm=250.0, filename=str(out))
        reloaded = load_normal_molecule_set(str(out))

        np.testing.assert_allclose(reloaded[:, :3], coords, rtol=1e-5)

    def test_stored_values_are_in_pixel_units(self, tmp_path):
        """The on-disk datatable holds pixels, not nanometres."""
        import h5py
        out = tmp_path / "units.h5"
        coords = np.array([[1000.0, 2000.0, 3000.0]])

        save_dataset_as_ms_h5(coords, np.array([0]), 100.0, pixelsize_z_nm=200.0,
                              filename=str(out))

        with h5py.File(str(out), "r") as f:
            row = f["molecule_set_data"]["datatable"][0]
            assert row["X_POS_PIXELS"] == pytest.approx(10.0)   # 1000 nm / 100
            assert row["Y_POS_PIXELS"] == pytest.approx(20.0)   # 2000 nm / 100
            assert row["Z_POS_PIXELS"] == pytest.approx(15.0)   # 3000 nm / 200
            assert f["molecule_set_data"]["xy_pixel_size_um"][()] == pytest.approx(0.1)
            assert f["molecule_set_data"]["z_pixel_size_um"][()] == pytest.approx(0.2)

    def test_2d_input_gives_zero_z(self, tmp_path):
        """(N, 2) input has no z; it must not come back with a spurious offset."""
        out = tmp_path / "two_d.h5"
        coords = np.array([[100.0, 500.0], [200.0, 600.0]])

        save_dataset_as_ms_h5(coords, np.array([0, 1]), 100.0, filename=str(out))
        reloaded = load_normal_molecule_set(str(out))

        np.testing.assert_allclose(reloaded[:, :2], coords, rtol=1e-5)
        assert np.all(reloaded[:, 2] == 0.0)


class TestMoleculeSetPixelSizeVariants:
    """The pixel-size key naming changed between molecule-set format versions."""

    @staticmethod
    def write_with_keys(path, keys):
        import h5py
        with h5py.File(str(path), "w") as f:
            g = f.create_group("molecule_set_data")
            for key, value_um in keys.items():
                g[key] = value_um
            dtype = np.dtype([("X_POS_PIXELS", np.float32), ("Y_POS_PIXELS", np.float32),
                              ("Z_POS_PIXELS", np.float32), ("FRAME_NUMBER", np.int32),
                              ("PHOTONS", np.float32)])
            table = np.zeros(1, dtype=dtype)
            table["X_POS_PIXELS"] = 10.0
            table["Y_POS_PIXELS"] = 20.0
            table["Z_POS_PIXELS"] = 5.0
            g.create_dataset("datatable", data=table)

    def test_legacy_single_general_pixel_size(self, tmp_path):
        out = tmp_path / "legacy.h5"
        self.write_with_keys(out, {"pixel_size_um": 0.1})

        locs = load_normal_molecule_set(str(out))

        # one general pixel size applies to every axis
        np.testing.assert_allclose(locs[0, :3], [1000.0, 2000.0, 500.0], rtol=1e-5)

    def test_legacy_alternative_spelling(self, tmp_path):
        out = tmp_path / "legacy2.h5"
        self.write_with_keys(out, {"pixelsize_um": 0.1})

        locs = load_normal_molecule_set(str(out))

        np.testing.assert_allclose(locs[0, :3], [1000.0, 2000.0, 500.0], rtol=1e-5)

    def test_separate_xy_and_z_keys(self, tmp_path):
        out = tmp_path / "modern.h5"
        self.write_with_keys(out, {"xy_pixel_size_um": 0.1, "z_pixel_size_um": 0.2})

        locs = load_normal_molecule_set(str(out))

        np.testing.assert_allclose(locs[0, :3], [1000.0, 2000.0, 1000.0], rtol=1e-5)

    def test_specific_keys_win_over_the_general_one(self, tmp_path):
        """A file carrying both must use the axis-specific z value."""
        out = tmp_path / "both.h5"
        self.write_with_keys(out, {"pixel_size_um": 0.1,
                                   "xy_pixel_size_um": 0.1,
                                   "z_pixel_size_um": 0.2})

        locs = load_normal_molecule_set(str(out))

        # z must be 5 px * 200 nm, not 5 px * 100 nm
        assert locs[0, 2] == pytest.approx(1000.0, rel=1e-5)

    def test_xy_key_without_z_falls_back_to_xy(self, tmp_path):
        out = tmp_path / "partial.h5"
        self.write_with_keys(out, {"xy_pixel_size_um": 0.1})

        locs = load_normal_molecule_set(str(out))

        np.testing.assert_allclose(locs[0, :3], [1000.0, 2000.0, 500.0], rtol=1e-5)

    def test_missing_pixel_size_reports_clearly(self, tmp_path):
        out = tmp_path / "nokeys.h5"
        self.write_with_keys(out, {})

        with pytest.raises(KeyError, match="No pixel size found"):
            load_normal_molecule_set(str(out))

    def test_photon_bandpass_filters_localizations(self, tmp_path, drifting_dataset):
        locs, _ = drifting_dataset
        out = tmp_path / "photons.h5"
        # half the localizations get a low photon count
        amplitudes = np.where(np.arange(len(locs)) % 2 == 0, 100.0, 5000.0)

        save_dataset_as_ms_h5(locs[:, :3], locs[:, 3], 100.0,
                              amplitudes=amplitudes, filename=str(out))

        everything = load_normal_molecule_set(str(out))
        bright_only = load_normal_molecule_set(str(out), photon_bandpass=(1000, 10000))

        assert len(bright_only) < len(everything)
        assert len(bright_only) == np.count_nonzero(amplitudes > 1000)

    def test_extra_dict_is_stored(self, tmp_path, drifting_dataset):
        import h5py
        locs, _ = drifting_dataset
        out = tmp_path / "extra.h5"

        save_dataset_as_ms_h5(locs[:, :3], locs[:, 3], 100.0, filename=str(out),
                              extra_dict={"note": 42})

        with h5py.File(str(out), "r") as f:
            assert f["extra_data"]["note"][()] == 42


class TestCorrectAndSave:
    def test_applies_drift_to_a_csv(self, tmp_path, thunderstorm_csv, drifting_dataset):
        locs, gt_drift = drifting_dataset
        n_frames = int(locs[:, 3].max()) + 1
        out = tmp_path / "corrected.csv"

        # a drift table shaped the way comet_run_kd returns it: (F, 4)
        drift_table = np.column_stack([gt_drift, np.arange(n_frames)])

        correct_and_save_thunderstorm_csv(drift_table, str(thunderstorm_csv), str(out))
        assert out.exists()

        corrected = load_thunderstorm_csv(str(out))
        original = load_thunderstorm_csv(str(thunderstorm_csv))
        frame = original[:, 3].astype(int)

        expected = original[:, :3] - gt_drift[frame]
        np.testing.assert_allclose(corrected[:, :3], expected, atol=1e-2)

    def test_short_drift_table_reports_clearly(self, tmp_path, thunderstorm_csv, drifting_dataset):
        """A drift table from a different dataset must not index out of bounds."""
        locs, gt_drift = drifting_dataset
        truncated = np.column_stack([gt_drift, np.arange(len(gt_drift))])[:5]

        with pytest.raises(ValueError, match="drift table only covers"):
            correct_and_save_thunderstorm_csv(
                truncated, str(thunderstorm_csv), str(tmp_path / "out.csv"))
