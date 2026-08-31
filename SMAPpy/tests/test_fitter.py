"""Fitter correctness on simulated data with known ground truth.

Simulating the Gaussian models is not circular: the simulation uses the
analytic integrated-Gaussian expression written independently of the fitter's
Levenberg-Marquardt code, so bias and precision here test the fitter itself.
(The spline model is checked against real data instead, since simulating from
the same coefficients it fits would only prove self-consistency.)
"""
import numpy as np
import pytest
from scipy.special import erf

from smappy.psf import GaussianPSF


def simulate(x, y, photons, background, sigma, roisize=13, seed=None,
             sigma_y=None):
    """Poisson ROIs of an integrated 2D Gaussian; x is the column index."""
    sigma_y = sigma if sigma_y is None else sigma_y
    i = np.arange(roisize)

    def prof(centre, s):
        return 0.5 * (erf((i - centre + 0.5) / (np.sqrt(2) * s))
                      - erf((i - centre - 0.5) / (np.sqrt(2) * s)))

    x, y, photons, background = np.broadcast_arrays(
        np.asarray(x, float), np.asarray(y, float),
        np.asarray(photons, float), np.asarray(background, float))
    n = x.size
    out = np.empty((n, roisize, roisize), np.float32)
    for k in range(n):
        px = prof(x.flat[k], sigma)
        py = prof(y.flat[k], sigma_y)
        out[k] = background.flat[k] + photons.flat[k] * np.outer(py, px)
    if seed is not None:
        out = np.random.default_rng(seed).poisson(out).astype(np.float32)
    return out


def test_noiseless_fit_recovers_the_truth():
    x, y = 6.3, 5.7
    rois = simulate(x, y, 2000.0, 10.0, 1.3)
    res = GaussianPSF(sigma=1.0).fit(rois, iterations=100)
    p = GaussianPSF(sigma=1.0).unpack(res)
    assert p["x_roi"][0] == pytest.approx(x, abs=1e-3)
    assert p["y_roi"][0] == pytest.approx(y, abs=1e-3)
    assert p["photons"][0] == pytest.approx(2000.0, rel=2e-3)
    assert p["background"][0] == pytest.approx(10.0, rel=2e-2)
    assert p["sigma_pix"][0] == pytest.approx(1.3, rel=2e-3)


def test_x_and_y_are_not_swapped():
    """An asymmetric spot must come back with x and y the right way round."""
    rois = simulate(8.0, 3.0, 3000.0, 5.0, 1.0, sigma_y=2.0)
    model = GaussianPSF(sigma=1.0, elliptical=True)
    p = model.unpack(model.fit(rois, iterations=100))
    assert p["x_roi"][0] == pytest.approx(8.0, abs=0.02)
    assert p["y_roi"][0] == pytest.approx(3.0, abs=0.02)
    assert p["sigma_x_pix"][0] == pytest.approx(1.0, rel=0.05)
    assert p["sigma_y_pix"][0] == pytest.approx(2.0, rel=0.05)


def test_unbiased_and_precision_matches_crlb():
    """With Poisson noise: no bias, and scatter consistent with the CRLB."""
    n = 3000
    rng = np.random.default_rng(1)
    x = rng.uniform(5.5, 6.5, n)
    y = rng.uniform(5.5, 6.5, n)
    rois = simulate(x, y, 1000.0, 10.0, 1.2, seed=2)
    model = GaussianPSF(sigma=1.2)
    p = model.unpack(model.fit(rois, iterations=100))

    dx, dy = p["x_roi"] - x, p["y_roi"] - y
    assert abs(np.mean(dx)) < 0.01 and abs(np.mean(dy)) < 0.01  # unbiased
    # measured scatter must agree with the predicted uncertainty
    assert np.std(dx) == pytest.approx(np.median(p["x_err_pix"]), rel=0.15)
    assert np.std(dy) == pytest.approx(np.median(p["y_err_pix"]), rel=0.15)
    assert np.mean(p["photons"]) == pytest.approx(1000.0, rel=0.02)


def test_precision_improves_with_photons():
    """More photons must localize better, and by the predicted amount.

    The scatter is compared against the fitter's own CRLB rather than against
    1/sqrt(N): with a background present the background term of the uncertainty
    scales as 1/N, so precision improves *faster* than sqrt(N) at low signal.
    """
    n = 800
    rng = np.random.default_rng(3)
    x = rng.uniform(5.5, 6.5, n)
    model = GaussianPSF(sigma=1.2)
    spread, predicted = [], []
    for photons in (500.0, 5000.0):
        rois = simulate(x, 6.0, photons, 10.0, 1.2, seed=int(photons))
        p = model.unpack(model.fit(rois, iterations=100))
        spread.append(np.std(p["x_roi"] - x))
        predicted.append(np.median(p["x_err_pix"]))

    assert spread[1] < spread[0]  # more photons localize better
    for measured, crlb in zip(spread, predicted):
        assert measured == pytest.approx(crlb, rel=0.2)


def test_threads_give_the_same_answer():
    rng = np.random.default_rng(4)
    rois = simulate(rng.uniform(5.5, 6.5, 500), 6.0, 1500.0, 8.0, 1.2, seed=5)
    model = GaussianPSF(sigma=1.2)
    one = model.fit(rois, n_threads=1)
    many = model.fit(rois, n_threads=8)
    assert np.array_equal(one.theta, many.theta)
    assert np.array_equal(one.logl, many.logl)
