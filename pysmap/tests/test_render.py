"""Geometry, accumulation and the display path, on positions with known pixels."""
import numpy as np
import pytest

from smapfit import lut as luts
from smapfit.render import (FieldOfView, SigmaSettings, normalize, render,
                            render_histogram, to_rgb)


def test_pixel_grid_and_edges():
    fov = FieldOfView.from_range((0, 100), (0, 50), 10.0)
    assert fov.shape == (5, 10)
    assert fov.extent == (0.0, 100.0, 50.0, 0.0)
    # a position sits in the pixel whose half-open interval contains it
    px, py = fov.to_pixels([0.0, 9.999, 10.0, 99.9], [0.0, 0.0, 0.0, 49.9])
    assert np.floor(px).tolist() == [0, 0, 1, 9]

    img = render_histogram([5.0, 15.0, 15.0], [5.0, 45.0, 45.0], fov)
    assert img.n_locs == 3
    assert img.weight[0, 0] == 1 and img.weight[4, 1] == 2
    assert img.weight.sum() == 3


def test_localizations_outside_are_dropped():
    fov = FieldOfView.from_range((0, 10), (0, 10), 1.0)
    img = render_histogram([-0.1, 5.0, 10.0], [5.0, 5.0, 5.0], fov)
    assert img.n_locs == 1 and img.weight.sum() == 1


def test_weights_are_summed():
    fov = FieldOfView.from_range((0, 4), (0, 4), 1.0)
    img = render_histogram([0.5, 0.5, 2.5], [0.5, 0.5, 0.5], fov, weights=[2.0, 3.0, 1.0])
    assert img.weight[0, 0] == 5.0 and img.weight[0, 2] == 1.0


def test_fit_keeps_pixels_square_and_centred():
    fov = FieldOfView.fit((0, 1000), (0, 500), 100, 100)
    assert fov.pixelsize == 10.0                      # the larger range decides
    assert fov.x0 == 0.0 and fov.y0 == -250.0         # the shorter one is centred
    assert 0.5 * (fov.y0 + fov.y1) == 250.0


def test_zoom_keeps_the_cursor_fixed():
    fov = FieldOfView.from_range((0, 100), (0, 100), 1.0)
    zoomed = fov.zoomed(0.5, cx=25.0, cy=75.0)
    assert zoomed.pixelsize == 0.5 and zoomed.shape == fov.shape
    # the anchor keeps its position on the canvas; the view shrinks around it
    assert np.allclose(zoomed.to_pixels(25.0, 75.0), fov.to_pixels(25.0, 75.0))
    assert (zoomed.x1 - zoomed.x0) == 0.5 * (fov.x1 - fov.x0)


def test_hue_and_sum_composites_agree_below_saturation():
    """The two field-colour composites differ only where a pixel saturates."""
    fov = FieldOfView.from_range((0, 4), (0, 4), 1.0)
    x = np.array([0.5, 0.5, 2.5])
    y = np.array([0.5, 0.5, 0.5])
    colors = luts.colors([0.0, 1.0, 0.5], "jet", 0.0, 1.0)
    img = render_histogram(x, y, fov, colors=colors)

    imax = float(img.weight.max())                     # nothing clips
    hue = to_rgb(img, imax=imax, color_mode="hue")
    total = to_rgb(img, imax=imax, color_mode="sum")
    assert np.allclose(hue, total, atol=1e-6)

    # the mixed pixel is the mean of the two colours it holds
    assert np.allclose(img.color[0, 0] / img.weight[0, 0], colors[:2].mean(axis=0))


def test_intensity_lut_is_applied_after_normalization():
    fov = FieldOfView.from_range((0, 2), (0, 1), 1.0)
    img = render_histogram([0.5], [0.5], fov)
    rgb = to_rgb(img, lut="gray", imax=1.0)
    assert np.allclose(rgb[0, 0], 1.0) and np.allclose(rgb[0, 1], 0.0)


def test_normalize_falls_back_to_the_maximum():
    """A very sparse image has a zero quantile; SMAP then uses the maximum."""
    image = np.zeros((100, 100), np.float32)
    image[0, 0] = 7.0
    norm, imax = normalize(image, contrast=0.3)  # a quantile that lands on zero
    assert imax == 7.0 and norm.max() == 1.0


def test_luts_are_well_formed():
    for name in luts.names():
        table = luts.get(name)
        assert table.shape == (256, 3) and table.min() >= 0 and table.max() <= 1
    assert np.allclose(luts.get("gray", invert=True), luts.get("gray_inverted"))
    # values below/above the range clamp to the first/last colour
    ends = luts.colors([-5.0, 5.0], "jet", 0.0, 1.0)
    assert np.allclose(ends, luts.get("jet")[[0, -1]])


# ------------------------------------------------------- Gaussian render modes
def _random_locs(n=200, seed=0):
    rng = np.random.default_rng(seed)
    return (rng.uniform(-20, 320, n).astype(np.float32),
            rng.uniform(-20, 220, n).astype(np.float32),
            rng.uniform(1.0, 20.0, n).astype(np.float32))


@pytest.mark.parametrize("mode", ["hist", "fixed", "per_loc"])
def test_extension_matches_the_numpy_reference(mode):
    x, y, sigma = _random_locs()
    fov = FieldOfView.from_range((0, 300), (0, 200), 5.0)
    sig = {"hist": None, "fixed": 12.0, "per_loc": sigma}[mode]
    colors = luts.colors(np.linspace(0, 1, x.size), "jet", 0.0, 1.0)

    fast = render(x, y, fov, sigma=sig, colors=colors)
    reference = render(x, y, fov, sigma=sig, colors=colors, use_extension=False)
    assert np.allclose(fast.weight, reference.weight, atol=1e-6)
    assert np.allclose(fast.color, reference.color, atol=1e-6)
    assert fast.n_locs == reference.n_locs


def test_every_localization_contributes_its_full_intensity():
    """The kernel is normalised by its own sum, so N is conserved exactly."""
    x, y, sigma = _random_locs(n=100, seed=1)
    fov = FieldOfView.around(x, y, 5.0, margin=200.0)  # nothing near the border
    weights = np.full(x.size, 3.0, np.float32)
    for sig in (None, 12.0, sigma):
        img = render(x, y, fov, sigma=sig, weights=weights)
        assert img.weight.sum() == pytest.approx(3.0 * x.size, rel=1e-5)


def test_the_histogram_is_the_zero_sigma_limit():
    x, y, _ = _random_locs(n=50, seed=2)
    fov = FieldOfView.from_range((0, 300), (0, 200), 5.0)
    hist = render(x, y, fov)
    tiny = render(x, y, fov, sigma=1e-3)
    assert np.allclose(hist.weight, tiny.weight, atol=1e-6)


def test_the_kernel_is_a_pixel_integrated_gaussian():
    """Its second moment must be sigma^2 + 1/12: the variance of the pixel itself.

    A point-sampled kernel (SMAP's template) would give sigma^2 instead, which
    is what makes it wrong for the small sigmas the renderer mostly uses.
    """
    fov = FieldOfView.from_range((0, 100), (0, 100), 1.0)
    for sigma in (0.8, 3.0):
        img = render(np.array([50.5]), np.array([50.5]), fov, sigma=sigma,
                     roi_sigma=8.0)
        rows, cols = np.mgrid[0:fov.ny, 0:fov.nx]
        centers = cols + 0.5
        mean = (img.weight * centers).sum()
        var = (img.weight * (centers - mean) ** 2).sum()
        assert mean == pytest.approx(50.5, abs=1e-3)
        assert var == pytest.approx(sigma ** 2 + 1 / 12, rel=1e-3)


def test_colors_do_not_change_the_weight_plane():
    x, y, sigma = _random_locs(n=100, seed=3)
    fov = FieldOfView.from_range((0, 300), (0, 200), 5.0)
    colors = luts.colors(sigma, "jet", sigma.min(), sigma.max())
    plain = render(x, y, fov, sigma=sigma)
    colored = render(x, y, fov, sigma=sigma, colors=colors)
    assert np.array_equal(plain.weight, colored.weight)
    # a white localization puts its whole weight into every channel
    white = render(x, y, fov, sigma=sigma, colors=np.ones((x.size, 3), np.float32))
    assert np.allclose(white.color, plain.weight[..., None], atol=1e-6)


def test_threads_do_not_change_the_result():
    x, y, sigma = _random_locs(n=500, seed=4)
    fov = FieldOfView.from_range((0, 300), (0, 200), 2.0)
    serial = render(x, y, fov, sigma=sigma, n_threads=1)
    threaded = render(x, y, fov, sigma=sigma, n_threads=8)
    assert np.array_equal(serial.weight, threaded.weight)
    assert serial.n_locs == threaded.n_locs


def test_sigma_settings_floor_and_cap():
    precision = np.array([0.5, 5.0, 5.0, 5.0, 500.0], np.float32)
    sigma = SigmaSettings(factor=2.0, min_sigma=4.0, min_sigma_pixels=0.7,
                          max_factor=10.0).apply(precision, pixelsize=10.0)
    assert sigma[0] == 7.0                  # floor: 0.7 pixels beats min_sigma
    assert sigma[1] == 10.0                 # 5 * factor
    assert sigma[4] == 100.0                # capped at 10x the median (10)


def test_dynamic_contrast_saturates_the_stated_fraction():
    """p is the handle SMAP uses: 10^-p of the pixels reach full scale."""
    rng = np.random.default_rng(0)
    image = rng.gamma(0.3, 10.0, (600, 600)).astype(np.float32)
    for p in (2.0, 3.0, 4.0):
        norm, imax = normalize(image, contrast=p)
        assert (image >= imax).mean() == pytest.approx(10.0 ** -p, rel=0.15)
        assert norm.max() == 1.0 and norm.min() >= 0.0
    # a lower p is a brighter image, at every pixel
    assert (normalize(image, contrast=2.0)[0] >= normalize(image, contrast=4.0)[0]).all()


def test_an_explicit_imax_overrides_the_contrast():
    image = np.linspace(0, 10, 100, dtype=np.float32).reshape(10, 10)
    norm, imax = normalize(image, imax=5.0, contrast=2.0)
    assert imax == 5.0 and norm.max() == 1.0 and (norm >= 1.0).sum() == 50


def test_the_two_composites_differ_only_where_a_pixel_saturates():
    """Red over cyan: SMAP's composite adds to white, the hue one stays grey.

    Below saturation the two are algebraically identical -- hue divides by the
    weight and multiplies it back in -- so this is the whole of the difference.
    """
    fov = FieldOfView.from_range((0, 3), (0, 1), 1.0)
    red, cyan = [1.0, 0.0, 0.0], [0.0, 1.0, 1.0]
    img = render_histogram([0.5, 1.5, 1.5, 2.5], [0.5] * 4, fov,
                           colors=np.array([red, red, cyan, cyan], np.float32))
    assert img.weight[0].tolist() == [1.0, 2.0, 1.0]

    below = [to_rgb(img, imax=2.0, color_mode=m) for m in ("hue", "sum")]
    assert np.allclose(below[0], below[1])

    hue = to_rgb(img, imax=1.0, color_mode="hue")
    total = to_rgb(img, imax=1.0, color_mode="sum")
    assert np.allclose(total[0, 1], [1, 1, 1])       # additive: white
    assert np.allclose(hue[0, 1], [0.5, 0.5, 0.5])   # hue: grey at full brightness
    for rgb in (hue, total):                          # the single colours agree
        assert np.allclose(rgb[0, 0], red) and np.allclose(rgb[0, 2], cyan)
