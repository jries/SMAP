"""Metadata merging and Micro-Manager parsing rules."""
import pytest

from smappy.metadata import CameraMetadata
from smappy.io.tiff import _parse_roi
from smappy.io.cameras_mat import _eval_expr


def test_roi_separator_is_not_a_minus_sign():
    assert _parse_roi("208-273-200-200") == (208, 273, 200, 200)
    assert _parse_roi("nonsense") is None
    assert _parse_roi(None) is None


def test_matlab_expressions():
    assert _eval_expr("str2double(X)", "6.7") == 6.7
    assert _eval_expr("~strcmp(X,'Normal')", "Multiplication Gain") is True
    assert _eval_expr("~strcmp(X,'Normal')", "Normal") is False
    assert _eval_expr("~strcmp(X,'Conventional')", "Electron Multiplying") is True
    assert _eval_expr("strcmp(X,'On')", "On") is True


def test_overrides_do_not_clear_unspecified_fields():
    """A config that only sets the pixel size must not switch EM off."""
    from_file = CameraMetadata(conversion=6.7, offset=400.0, em_on=True, emgain=100)
    merged = from_file.merged_with(CameraMetadata(pixelsize_um=0.127))
    assert merged.em_on is True and merged.emgain == 100
    assert merged.pixelsize_um == 0.127
    assert merged.adu_to_photons == pytest.approx(0.067)
    assert merged.excess_noise == 2.0


def test_overrides_win_when_set():
    from_file = CameraMetadata(conversion=6.7, offset=400.0, em_on=True, emgain=100)
    merged = from_file.merged_with(CameraMetadata(offset=398.6, em_on=False))
    assert merged.offset == 398.6
    assert merged.is_em is False
    assert merged.adu_to_photons == 6.7
    assert merged.excess_noise == 1.0


def test_missing_values_are_reported():
    with pytest.raises(ValueError, match="pixelsize_um"):
        CameraMetadata(conversion=6.7, offset=400.0).require()
