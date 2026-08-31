"""Unit conversion, and the HDF5 round trip."""
import numpy as np
import pytest

from smappy.locs import Localizations, to_nm
from smappy.io.hdf5 import LocalizationWriter, load_localizations


def _table(n=100, seed=0):
    rng = np.random.default_rng(seed)
    return Localizations({
        "frame": np.arange(n, dtype=np.int64),
        "x_pix": rng.uniform(0, 200, n).astype(np.float32),
        "y_pix": rng.uniform(0, 200, n).astype(np.float32),
        "x_err_pix": rng.uniform(0.01, 0.2, n).astype(np.float32),
        "y_err_pix": rng.uniform(0.01, 0.2, n).astype(np.float32),
        "photons": rng.uniform(500, 5000, n).astype(np.float32),
        "z_nm": rng.uniform(-400, 400, n).astype(np.float32),
    }, {"units": "pixel"})


def test_pixels_are_primary_and_nm_is_derived():
    locs = _table()
    nm = to_nm(locs, 127.0)
    assert "x_nm" in nm and "x_pix" not in nm
    assert np.allclose(nm["x_nm"], locs["x_pix"] * 127.0)
    assert np.allclose(nm["x_err_nm"], locs["x_err_pix"] * 127.0)
    # z is already nm and must not be scaled again
    assert np.array_equal(nm["z_nm"], locs["z_nm"])
    # photons are not a length
    assert np.array_equal(nm["photons"], locs["photons"])
    assert nm.metadata["units"] == "nm" and nm.metadata["pixelsize_nm"] == 127.0


def test_keep_pixels_gives_both():
    both = to_nm(_table(), 100.0, keep_pixels=True)
    assert "x_pix" in both and "x_nm" in both
    assert both.metadata["units"] == "pixel+nm"


def test_converting_twice_is_a_no_op():
    once = to_nm(_table(), 127.0)
    assert to_nm(once, 127.0) is once


def test_hdf5_round_trip(tmp_path):
    locs = _table()
    path = tmp_path / "locs.h5"
    with LocalizationWriter(path, {"source": "test"}) as writer:
        writer.append(locs)
    back = load_localizations(path)
    assert len(back) == len(locs)
    for name in locs.keys():
        assert np.allclose(back[name], locs[name])
    assert back.metadata["source"] == "test"
    assert back.metadata["units"] == "pixel"  # taken from the table itself


def test_streaming_equals_one_shot(tmp_path):
    """Appending in blocks must give the same file as writing in one go."""
    locs = _table(300, seed=1)
    stream, single = tmp_path / "a.h5", tmp_path / "b.h5"
    with LocalizationWriter(stream) as writer:
        for start in range(0, 300, 37):
            writer.append(locs[slice(start, start + 37)])
    with LocalizationWriter(single) as writer:
        writer.append(locs)
    a, b = load_localizations(stream), load_localizations(single)
    assert len(a) == len(b) == 300
    for name in locs.keys():
        assert np.array_equal(a[name], b[name])


def test_partial_file_is_readable(tmp_path):
    """Results must be on disk before the writer is closed."""
    path = tmp_path / "partial.h5"
    writer = LocalizationWriter(path)
    writer.append(_table(50))
    assert len(load_localizations(path)) == 50  # readable while still open
    writer.append(_table(50, seed=2))
    writer.close()
    assert len(load_localizations(path)) == 100


def test_mismatched_columns_are_refused(tmp_path):
    writer = LocalizationWriter(tmp_path / "x.h5")
    writer.append(_table(10))
    other = _table(10)
    del other.columns["photons"]
    with pytest.raises(ValueError, match="columns changed"):
        writer.append(other)
    writer.close()


def test_prefetch_preserves_order_and_content():
    from smappy.pipeline import prefetch
    source = [(i, np.full((2, 3), i, np.float32)) for i in range(20)]
    out = list(prefetch(iter(source), depth=3))
    assert [i for i, _ in out] == list(range(20))
    assert all(np.array_equal(a, b) for (_, a), (_, b) in zip(out, source))


def test_prefetch_reraises_reader_errors():
    from smappy.pipeline import prefetch

    def broken():
        yield (0, np.zeros((2, 2), np.float32))
        raise RuntimeError("disk went away")

    with pytest.raises(RuntimeError, match="disk went away"):
        list(prefetch(broken(), depth=2))
