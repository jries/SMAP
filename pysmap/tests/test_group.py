"""Linking and combining, on tables whose groups are known by construction."""
import numpy as np
import pytest

from smapfit.group import COMBINE_MODES, GroupSettings, combine, connect, group
from smapfit.locs import Localizations


def _runs(spec, dx_jitter=0.0, seed=0):
    """A table from (x, y, first_frame, length) runs; one emitter per run."""
    rng = np.random.default_rng(seed)
    x, y, frame, truth = [], [], [], []
    for i, (px, py, f0, n) in enumerate(spec):
        for k in range(n):
            x.append(px + rng.normal(0, dx_jitter))
            y.append(py + rng.normal(0, dx_jitter))
            frame.append(f0 + k)
            truth.append(i)
    order = rng.permutation(len(x))          # input order must not matter
    return (np.array(x)[order], np.array(y)[order],
            np.array(frame)[order], np.array(truth)[order])


def _same_partition(a, b):
    """Whether two labellings define the same grouping, ignoring the labels."""
    return len(set(zip(np.asarray(a).tolist(), np.asarray(b).tolist()))) == \
        len(np.unique(a)) == len(np.unique(b))


def test_consecutive_frames_are_linked():
    x, y, frame, truth = _runs([(100, 100, 0, 5), (300, 300, 0, 3),
                                (100, 100, 20, 4)], dx_jitter=2.0)
    ids = connect(x, y, frame, dx=20.0, dt=0)
    assert _same_partition(ids, truth)
    assert ids.min() == 1


def test_a_gap_larger_than_dt_starts_a_new_group():
    # one emitter, on for 3 frames, dark for 2, on for 3
    x = np.full(6, 100.0)
    y = np.full(6, 100.0)
    frame = np.array([0, 1, 2, 5, 6, 7])
    assert len(np.unique(connect(x, y, frame, dx=20.0, dt=1))) == 2
    assert len(np.unique(connect(x, y, frame, dx=20.0, dt=2))) == 1


def test_the_search_region_is_a_box():
    """SMAP tests |dx| and |dy| separately, so the corner of the box links."""
    frame = np.array([0, 1])
    inside = connect(np.array([0.0, 9.0]), np.array([0.0, 9.0]), frame, 10.0, 0)
    outside = connect(np.array([0.0, 11.0]), np.array([0.0, 0.0]), frame, 10.0, 0)
    assert len(np.unique(inside)) == 1      # 12.7 away, but inside the box
    assert len(np.unique(outside)) == 2


def test_linking_does_not_cross_blocks():
    x = np.full(4, 100.0)
    y = np.full(4, 100.0)
    frame = np.array([0, 1, 2, 3])
    assert len(np.unique(connect(x, y, frame, 20.0, 1))) == 1
    blocks = np.array([0, 0, 1, 1])
    ids = connect(x, y, frame, 20.0, 1, blocks=blocks)
    assert len(np.unique(ids)) == 2 and sorted(np.unique(ids)) == [1, 2]


def test_every_localization_gets_exactly_one_group():
    """The bug fixed in the C port: SMAP left the first and last unassigned."""
    x, y, frame, _ = _runs([(i * 50.0, 0.0, i % 7, 1 + i % 4) for i in range(60)],
                           dx_jitter=1.0, seed=3)
    ids = connect(x, y, frame, dx=20.0, dt=1)
    assert (ids > 0).all()
    assert np.array_equal(np.unique(ids), np.arange(1, ids.max() + 1))


def _table(n=12):
    rng = np.random.default_rng(1)
    return Localizations({
        "x_nm": np.repeat([100.0, 500.0, 900.0], 4).astype(np.float32),
        "y_nm": np.zeros(n, np.float32),
        "loc_precision_nm": rng.uniform(5, 20, n).astype(np.float32),
        "photons": rng.uniform(1000, 5000, n).astype(np.float32),
        "photons_err": rng.uniform(50, 200, n).astype(np.float32),
        "logl_rel": rng.normal(-1, 0.3, n).astype(np.float32),
        "frame": np.arange(n, dtype=np.int64) % 4,
    }, {"units": "nm"})


def test_combine_follows_the_rules_per_column():
    locs = _table()
    gi = np.repeat([1, 2, 3], 4)
    grouped = combine(locs, gi)
    assert len(grouped) == 3
    assert np.array_equal(grouped["n_in_group"], [4, 4, 4])

    w = 1.0 / np.asarray(locs["loc_precision_nm"], np.float64) ** 2
    for g in range(3):
        s = slice(4 * g, 4 * g + 4)
        assert grouped["x_nm"][g] == pytest.approx(
            np.average(locs["x_nm"][s], weights=w[s]), rel=1e-5)
        assert grouped["photons"][g] == pytest.approx(
            np.sum(locs["photons"][s]), rel=1e-5)
        assert grouped["loc_precision_nm"][g] == pytest.approx(
            1 / np.sqrt(np.sum(1 / locs["loc_precision_nm"][s] ** 2)), rel=1e-5)
        # the error of a *sum* adds in quadrature, not by the precision rule
        assert grouped["photons_err"][g] == pytest.approx(
            np.sqrt(np.sum(locs["photons_err"][s] ** 2)), rel=1e-5)
        assert grouped["logl_rel"][g] == pytest.approx(max(locs["logl_rel"][s]))
        assert grouped["frame"][g] == min(locs["frame"][s])


def test_grouping_improves_precision_and_conserves_photons():
    locs = _table()
    gi = np.repeat([1, 2, 3], 4)
    grouped = combine(locs, gi)
    assert np.asarray(grouped["photons"]).sum() == pytest.approx(
        np.asarray(locs["photons"], np.float64).sum(), rel=1e-5)
    for g in range(3):
        best = locs["loc_precision_nm"][4 * g:4 * g + 4].min()
        assert grouped["loc_precision_nm"][g] < best


def test_the_meaningless_columns_are_dropped_not_guessed():
    """A group's raw logl is a sum over its members; no reduction is right."""
    locs = _table()
    locs.columns["logl"] = np.full(len(locs), -100.0, np.float32)
    grouped = combine(locs, np.repeat([1, 2, 3], 4))
    assert "logl" not in grouped and "logl_rel" in grouped
    assert COMBINE_MODES["logl_rel"] == "max"


def test_group_end_to_end_and_metadata():
    locs = _table()
    grouped, gi = group(locs, GroupSettings(dx=50.0, dt=1))
    assert len(grouped) == len(np.unique(gi))
    assert np.asarray(grouped["n_in_group"]).sum() == len(locs)
    assert grouped.metadata["grouped"] is True
    assert grouped.metadata["units"] == "nm"


def test_summed_errors_stay_consistent_with_shot_noise():
    """std(N) = sqrt(N) must survive grouping: std(N1+N2) = sqrt(N1+N2).

    This is what picks quadrature over SMAP's precision rule for the error of a
    summed quantity -- the identity holds exactly, not approximately.
    """
    photons = np.array([1000.0, 2000.0, 4000.0, 300.0], np.float32)
    locs = Localizations({
        "x_nm": np.zeros(4, np.float32), "y_nm": np.zeros(4, np.float32),
        "photons": photons, "photons_err": np.sqrt(photons).astype(np.float32),
        "loc_precision_nm": np.full(4, 10.0, np.float32),
        "frame": np.arange(4, dtype=np.int64),
    }, {"units": "nm"})
    grouped = combine(locs, np.ones(4, np.int64))
    assert grouped["photons"][0] == pytest.approx(photons.sum(), rel=1e-5)
    assert grouped["photons_err"][0] == pytest.approx(np.sqrt(photons.sum()),
                                                      rel=1e-5)
    # SMAP's rule would give this instead, six times too small
    smap = 1 / np.sqrt(np.sum(1 / np.sqrt(photons.astype(np.float64)) ** 2))
    assert smap < 0.2 * grouped["photons_err"][0]


def test_each_coordinate_is_weighted_by_its_own_error():
    """Under astigmatism x_err and y_err diverge, so the pooled weight is wrong."""
    locs = Localizations({
        "x_nm": np.array([-100.0, 0.0, 40.0], np.float32),
        "y_nm": np.array([40.0, 0.0, -100.0], np.float32),
        "x_err_nm": np.array([50.0, 5.0, 50.0], np.float32),
        "y_err_nm": np.array([5.0, 50.0, 5.0], np.float32),
        "loc_precision_nm": np.array([35.4, 35.4, 35.4], np.float32),
        "frame": np.arange(3, dtype=np.int64),
    }, {"units": "nm"})
    grouped = combine(locs, np.ones(3, np.int64))

    for name, err in (("x_nm", "x_err_nm"), ("y_nm", "y_err_nm")):
        values = np.asarray(locs[name], np.float64)
        own = np.average(values, weights=1 / np.asarray(locs[err], np.float64) ** 2)
        pooled = values.mean()          # the pooled precision is equal here
        assert grouped[name][0] == pytest.approx(own, abs=1e-3)
        assert abs(own - pooled) > 1.0  # the two really do differ


def test_z_is_weighted_by_its_own_error_when_the_table_has_one():
    z = np.array([-100.0, 0.0, 40.0], np.float32)   # asymmetric, so the
    # two weightings cannot coincide by symmetry
    z_err = np.array([50.0, 5.0, 50.0], np.float32)      # the middle one wins
    locs = Localizations({
        "x_nm": np.zeros(3, np.float32), "y_nm": np.zeros(3, np.float32),
        "z_nm": z, "z_err_nm": z_err,
        "loc_precision_nm": np.array([5.0, 50.0, 5.0], np.float32),
        "frame": np.arange(3, dtype=np.int64),
    }, {"units": "nm"})
    grouped = combine(locs, np.ones(3, np.int64))

    by_z_err = np.average(z.astype(np.float64), weights=1 / z_err.astype(np.float64) ** 2)
    by_lateral = np.average(z.astype(np.float64),
                            weights=1 / locs["loc_precision_nm"].astype(np.float64) ** 2)
    assert grouped["z_nm"][0] == pytest.approx(by_z_err, abs=1e-3)
    assert abs(by_z_err - by_lateral) > 1.0        # the two really do differ here
    assert grouped["z_err_nm"][0] == pytest.approx(
        1 / np.sqrt(np.sum(1 / z_err.astype(np.float64) ** 2)), rel=1e-5)
