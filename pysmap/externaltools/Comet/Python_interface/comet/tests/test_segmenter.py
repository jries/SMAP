"""Tests for temporal segmentation (all three modes plus edge cases)."""

import numpy as np
import pytest

from comet.core.segmenter import (
    SegmentationResult,
    segment_by_frame_windows,
    segment_by_num_locs_per_window,
    segment_by_num_windows,
    segmentation_wrapper,
)


def frames(n_frames, locs_per_frame):
    return np.repeat(np.arange(n_frames), locs_per_frame)


class TestFrameWindows:
    """Mode 2: fixed number of frames per window."""

    def test_segment_count_and_assignment(self):
        loc_frames = frames(n_frames=20, locs_per_frame=5)
        result = segment_by_frame_windows(loc_frames, n_frames_per_window=4)

        assert result.n_segments == 5
        assert set(np.unique(result.loc_segments)) == set(range(5))
        # every localization is kept when no downsampling cap is given
        assert result.loc_valid.all()

    def test_frames_map_to_expected_window(self):
        loc_frames = frames(n_frames=10, locs_per_frame=2)
        result = segment_by_frame_windows(loc_frames, n_frames_per_window=5)

        # frames 0-4 -> segment 0, frames 5-9 -> segment 1
        assert np.all(result.loc_segments[loc_frames < 5] == 0)
        assert np.all(result.loc_segments[loc_frames >= 5] == 1)

    def test_partial_trailing_window(self):
        # 10 frames in windows of 3 -> 4 segments, the last holding one frame
        result = segment_by_frame_windows(frames(10, 4), n_frames_per_window=3)
        assert result.n_segments == 4

    def test_center_frames_are_within_their_window(self):
        loc_frames = frames(n_frames=30, locs_per_frame=3)
        result = segment_by_frame_windows(loc_frames, n_frames_per_window=6)

        for seg in range(result.n_segments):
            member_frames = loc_frames[result.loc_segments == seg]
            assert member_frames.min() <= result.center_frames[seg] <= member_frames.max()

    def test_window_larger_than_dataset_yields_single_segment(self):
        result = segment_by_frame_windows(frames(5, 2), n_frames_per_window=1000)
        assert result.n_segments == 1
        assert np.all(result.loc_segments == 0)


class TestLocsPerWindow:
    """Mode 1: minimum number of localizations per window."""

    def test_each_segment_meets_the_minimum(self):
        loc_frames = frames(n_frames=40, locs_per_frame=10)  # 400 locs
        result = segment_by_num_locs_per_window(loc_frames, min_n_locs_per_window=50,
                                                return_param_dict=True)

        counts = result.out_dict["locs_per_segment"]
        assert result.n_segments > 1
        # the trailing segment absorbs the remainder, so it may exceed but not undershoot
        assert counts.min() >= 50

    def test_all_locs_assigned_to_a_segment(self):
        loc_frames = frames(n_frames=30, locs_per_frame=7)
        result = segment_by_num_locs_per_window(loc_frames, min_n_locs_per_window=35)
        # -1 is the "unassigned" sentinel
        assert not np.any(result.loc_segments == -1)


class TestNumWindows:
    """Mode 0: requested number of windows."""

    @pytest.mark.parametrize("n_windows", [2, 4, 8])
    def test_approximately_the_requested_window_count(self, n_windows):
        loc_frames = frames(n_frames=64, locs_per_frame=8)
        result = segment_by_num_windows(loc_frames, n_windows=n_windows)
        # derived via a locs-per-window threshold, so allow one segment of slack
        assert abs(result.n_segments - n_windows) <= 1


class TestWrapperDispatch:
    def test_mode_0_matches_num_windows(self):
        loc_frames = frames(20, 5)
        assert (segmentation_wrapper(loc_frames, 4, segmentation_mode=0).n_segments
                == segment_by_num_windows(loc_frames, 4).n_segments)

    def test_mode_1_matches_locs_per_window(self):
        loc_frames = frames(20, 5)
        assert (segmentation_wrapper(loc_frames, 25, segmentation_mode=1).n_segments
                == segment_by_num_locs_per_window(loc_frames, 25).n_segments)

    def test_mode_2_matches_frame_windows(self):
        loc_frames = frames(20, 5)
        assert (segmentation_wrapper(loc_frames, 4, segmentation_mode=2).n_segments
                == segment_by_frame_windows(loc_frames, 4).n_segments)

    def test_mode_2_is_the_default(self):
        loc_frames = frames(20, 5)
        assert (segmentation_wrapper(loc_frames, 4).n_segments
                == segment_by_frame_windows(loc_frames, 4).n_segments)

    def test_returns_segmentation_result(self):
        result = segmentation_wrapper(frames(10, 2), 5)
        assert isinstance(result, SegmentationResult)
        assert result.out_dict is None  # only populated on request


class TestDownsampling:
    def test_max_locs_per_segment_caps_valid_locs(self):
        loc_frames = frames(n_frames=20, locs_per_frame=50)  # 250 locs per 5-frame window
        result = segment_by_frame_windows(loc_frames, n_frames_per_window=5,
                                         max_locs_per_segment=40, return_param_dict=True)

        for seg in range(result.n_segments):
            kept = np.count_nonzero((result.loc_segments == seg) & result.loc_valid)
            assert kept <= 40
        assert result.out_dict["n_locs_valid"] < result.out_dict["n_locs"]

    def test_cap_above_segment_size_keeps_everything(self):
        result = segment_by_frame_windows(frames(10, 5), n_frames_per_window=5,
                                         max_locs_per_segment=10 ** 6)
        assert result.loc_valid.all()


class TestEdgeCases:
    def test_single_frame(self):
        result = segment_by_frame_windows(np.zeros(20, dtype=int), n_frames_per_window=4)
        assert result.n_segments == 1
        assert result.loc_valid.all()

    def test_sparse_and_unsorted_frames(self):
        # gaps between frames and arbitrary ordering must both be tolerated
        loc_frames = np.array([9, 0, 5, 9, 0, 21, 5, 21, 13, 13])
        result = segment_by_frame_windows(loc_frames, n_frames_per_window=2)

        assert not np.any(result.loc_segments == -1)
        # identical frames must land in the same segment regardless of position
        for frame in np.unique(loc_frames):
            segs = result.loc_segments[loc_frames == frame]
            assert len(np.unique(segs)) == 1

    def test_param_dict_totals_are_consistent(self):
        loc_frames = frames(24, 6)
        result = segment_by_frame_windows(loc_frames, 4, return_param_dict=True)
        d = result.out_dict

        assert d["n_locs"] == len(loc_frames)
        assert d["n_locs_valid"] + d["n_locs_invalid"] == d["n_locs"]
        assert d["n_segments"] == result.n_segments
        assert d["frames_per_window"] == 4
        assert len(d["start_frames"]) == len(d["end_frames"]) == result.n_segments
