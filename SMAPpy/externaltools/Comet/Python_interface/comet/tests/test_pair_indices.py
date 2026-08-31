"""Tests for the KD-tree neighbour-pair search.

The pair list drives the whole cost function, so it is checked against a
brute-force reference on inputs small enough to enumerate exhaustively.
"""

import numpy as np
import pytest

from comet.core.pair_indices import pair_indices_kdtree


def brute_force_pairs(coords, distance):
    """All (i, j) with i < j and |coords[i] - coords[j]| <= distance."""
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    i, j = np.triu_indices(len(coords), k=1)
    keep = dist[i, j] <= distance
    return set(zip(i[keep].tolist(), j[keep].tolist()))


def as_pair_set(idx_i, idx_j):
    return set(zip(np.asarray(idx_i).tolist(), np.asarray(idx_j).tolist()))


class TestCorrectness:
    @pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
    def test_matches_brute_force_3d(self, seed):
        rng = np.random.default_rng(seed)
        coords = rng.random((60, 3)) * 100.0
        distance = 25.0

        idx_i, idx_j, ok = pair_indices_kdtree(coords, distance)

        assert ok
        assert as_pair_set(idx_i, idx_j) == brute_force_pairs(coords, distance)

    @pytest.mark.parametrize("distance", [5.0, 20.0, 50.0, 200.0])
    def test_matches_brute_force_across_radii(self, distance):
        rng = np.random.default_rng(7)
        coords = rng.random((50, 3)) * 100.0

        idx_i, idx_j, ok = pair_indices_kdtree(coords, distance)

        assert ok
        assert as_pair_set(idx_i, idx_j) == brute_force_pairs(coords, distance)

    def test_known_geometry(self):
        # points on a line 10 nm apart: radius 15 links only adjacent neighbours
        coords = np.array([[0.0, 0, 0], [10.0, 0, 0], [20.0, 0, 0], [30.0, 0, 0]])
        idx_i, idx_j, ok = pair_indices_kdtree(coords, 15.0)

        assert ok
        assert as_pair_set(idx_i, idx_j) == {(0, 1), (1, 2), (2, 3)}

    def test_radius_is_inclusive(self):
        coords = np.array([[0.0, 0, 0], [10.0, 0, 0]])
        _, _, ok = pair_indices_kdtree(coords, 10.0)
        idx_i, idx_j, ok = pair_indices_kdtree(coords, 10.0)

        assert ok
        assert as_pair_set(idx_i, idx_j) == {(0, 1)}


class TestConventions:
    def test_indices_are_ordered_and_never_self_paired(self):
        rng = np.random.default_rng(11)
        coords = rng.random((80, 3)) * 50.0

        idx_i, idx_j, ok = pair_indices_kdtree(coords, 20.0)

        assert ok
        assert np.all(np.asarray(idx_i) < np.asarray(idx_j))

    def test_no_duplicate_pairs(self):
        rng = np.random.default_rng(12)
        coords = rng.random((80, 3)) * 50.0

        idx_i, idx_j, ok = pair_indices_kdtree(coords, 20.0)

        assert ok
        pairs = as_pair_set(idx_i, idx_j)
        assert len(pairs) == len(np.asarray(idx_i))

    def test_dtype_is_int32_and_contiguous(self):
        # the CUDA kernels index device arrays with int32
        rng = np.random.default_rng(13)
        coords = rng.random((40, 3)) * 50.0

        idx_i, idx_j, ok = pair_indices_kdtree(coords, 20.0)

        assert ok
        for arr in (idx_i, idx_j):
            assert arr.dtype == np.int32
            assert arr.flags["C_CONTIGUOUS"]

    def test_does_not_mutate_input(self):
        rng = np.random.default_rng(14)
        coords = rng.random((40, 3)) * 50.0
        original = coords.copy()

        pair_indices_kdtree(coords, 20.0)

        np.testing.assert_array_equal(coords, original)

    def test_estimate_pairs_does_not_mutate_input(self):
        """estimate_pairs shifts coordinates to the origin; on a caller's array
        that silently moved their data."""
        from comet.core.pair_indices import estimate_pairs

        rng = np.random.default_rng(16)
        coords = rng.random((60, 3)) * 500.0 + 1000.0
        original = coords.copy()

        estimate_pairs(coords, 50.0)

        np.testing.assert_array_equal(coords, original)

    def test_lex_floor_does_not_mutate_input(self):
        from comet.core.pair_indices import pair_indices_lex_floor_asymmetric

        rng = np.random.default_rng(17)
        coords = rng.random((60, 3)) * 500.0 + 1000.0
        original = coords.copy()

        pair_indices_lex_floor_asymmetric(coords, 50.0)

        np.testing.assert_array_equal(coords, original)


class TestEdgeCases:
    def test_no_pairs_within_radius(self):
        coords = np.array([[0.0, 0, 0], [1000.0, 0, 0], [2000.0, 0, 0]])
        idx_i, idx_j, ok = pair_indices_kdtree(coords, 10.0)

        assert ok
        assert len(idx_i) == 0 and len(idx_j) == 0

    def test_coincident_points_pair_up(self):
        coords = np.zeros((4, 3))
        idx_i, idx_j, ok = pair_indices_kdtree(coords, 1.0)

        assert ok
        # every unordered pair of the four identical points
        assert len(idx_i) == 6

    def test_two_dimensional_input(self):
        rng = np.random.default_rng(15)
        coords = rng.random((40, 2)) * 100.0
        distance = 25.0

        idx_i, idx_j, ok = pair_indices_kdtree(coords, distance)

        assert ok
        assert as_pair_set(idx_i, idx_j) == brute_force_pairs(coords, distance)
