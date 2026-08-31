"""Tests for per-segment -> per-frame drift interpolation."""

import numpy as np
import pytest

from comet.core.interpolation import interpolate_drift


def drift_from(fn, center_frames):
    """Build an (M, 3) drift array by evaluating fn on the segment centers."""
    values = fn(center_frames)
    return np.column_stack([values, 0.5 * values, -0.25 * values])


class TestExactness:
    """A well-behaved interpolant must reproduce its input at the knots."""

    @pytest.mark.parametrize("method", ["cubic", "linear"])
    def test_passes_through_control_points(self, method):
        centers = np.arange(0, 100, 10, dtype=float)
        drift = drift_from(lambda t: 3.0 * np.sin(t / 15.0), centers)

        out = interpolate_drift(centers, drift, centers, method=method)
        np.testing.assert_allclose(out, drift, atol=1e-9)

    def test_cubic_reproduces_a_cubic_polynomial(self):
        # a cubic spline is exact for data drawn from a cubic
        centers = np.linspace(0, 90, 10)
        query = np.arange(0, 91, dtype=float)
        poly = lambda t: 1e-4 * t ** 3 - 0.02 * t ** 2 + 0.5 * t + 3.0

        out = interpolate_drift(centers, drift_from(poly, centers), query, method="cubic")
        np.testing.assert_allclose(out[:, 0], poly(query), atol=1e-6)

    def test_linear_reproduces_a_straight_line(self):
        centers = np.linspace(0, 50, 6)
        query = np.arange(0, 51, dtype=float)
        line = lambda t: 2.5 * t - 7.0

        out = interpolate_drift(centers, drift_from(line, centers), query, method="linear")
        np.testing.assert_allclose(out[:, 0], line(query), atol=1e-9)


class TestShapeAndAxes:
    @pytest.mark.parametrize("method", ["cubic", "linear", "catmull-rom"])
    def test_output_shape(self, method):
        centers = np.linspace(0, 100, 12)
        query = np.arange(0, 101, dtype=float)
        out = interpolate_drift(centers, drift_from(np.sin, centers), query, method=method)
        assert out.shape == (len(query), 3)

    @pytest.mark.parametrize("method", ["cubic", "linear"])
    def test_axes_are_independent(self, method):
        # y is 0.5x and z is -0.25x by construction; interpolation must preserve that
        centers = np.linspace(0, 80, 9)
        query = np.arange(0, 81, dtype=float)
        out = interpolate_drift(centers, drift_from(lambda t: np.sin(t / 9.0), centers),
                                query, method=method)

        np.testing.assert_allclose(out[:, 1], 0.5 * out[:, 0], atol=1e-9)
        np.testing.assert_allclose(out[:, 2], -0.25 * out[:, 0], atol=1e-9)

    @pytest.mark.parametrize("method", ["cubic", "linear", "catmull-rom"])
    def test_no_nans_inside_the_knot_range(self, method):
        centers = np.linspace(0, 100, 11)
        # catmull-rom only defines the interior spans, so stay well inside
        query = np.arange(20, 81, dtype=float)
        out = interpolate_drift(centers, drift_from(np.cos, centers), query, method=method)
        assert np.isfinite(out).all()


class TestCatmullRom:
    def test_result_is_float_for_an_integer_frame_range(self):
        """comet_run_kd passes an int frame range; the output must not be truncated.

        np.zeros_like(x_interp) inherited int64 there, rounding every drift value
        to whole nanometres.
        """
        centers = np.linspace(0, 100, 11)
        drift = drift_from(lambda t: 3.7 * np.sin(t / 15.0), centers)
        frames_int = np.arange(0, 101, dtype=int)

        out = interpolate_drift(centers, drift, frames_int, method="catmull-rom")

        assert out.dtype.kind == "f"
        # sub-nanometre structure must survive
        assert not np.allclose(out[:, 0], np.round(out[:, 0]))

    def test_edges_are_clamped_not_zeroed(self):
        """Frames outside the interior spans must not read as 'no drift'."""
        centers = np.linspace(0, 100, 11)
        offset = 50.0                       # far from zero, so zeros are obvious
        drift = drift_from(lambda t: offset + 0.1 * t, centers)
        frames_int = np.arange(0, 101, dtype=int)

        out = interpolate_drift(centers, drift, frames_int, method="catmull-rom")

        assert np.isfinite(out).all()
        # nothing should have collapsed to zero at the ends
        assert np.abs(out[:, 0]).min() > offset / 2
        # the clamped regions hold the nearest in-range knot value
        np.testing.assert_allclose(out[: int(centers[1]), 0], drift[1, 0])
        np.testing.assert_allclose(out[int(centers[-2]) + 1:, 0], drift[-2, 0])

    def test_matches_a_straight_line_on_interior_spans(self):
        centers = np.linspace(0, 100, 11)
        query = np.arange(10, 91, dtype=float)
        line = lambda t: 0.3 * t + 1.0

        out = interpolate_drift(centers, drift_from(line, centers), query, method="catmull-rom")
        np.testing.assert_allclose(out[:, 0], line(query), atol=1e-6)

    def test_requires_at_least_four_control_points(self):
        centers = np.array([0.0, 10.0, 20.0])
        with pytest.raises(ValueError, match="at least 4 points"):
            interpolate_drift(centers, drift_from(np.sin, centers),
                              np.arange(0, 21, dtype=float), method="catmull-rom")


class TestErrors:
    def test_unknown_method_raises(self):
        centers = np.linspace(0, 50, 6)
        with pytest.raises(ValueError, match="Unknown interpolation method"):
            interpolate_drift(centers, drift_from(np.sin, centers),
                              np.arange(0, 51, dtype=float), method="nope")
