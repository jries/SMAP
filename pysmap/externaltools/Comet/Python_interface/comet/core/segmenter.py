import numpy as np
from dataclasses import dataclass
from typing import Optional, Dict


@dataclass
class SegmentationResult:
    """Container for segmentation results."""
    loc_segments: np.ndarray
    loc_valid: np.ndarray
    center_frames: np.ndarray
    n_segments: int
    out_dict: Optional[Dict] = None


def _group_by_frame(loc_frames: np.ndarray):
    """Returns a dict {frame_number: indices_in_loc_frames} efficiently."""
    sort_idx = np.argsort(loc_frames)
    sorted_frames = loc_frames[sort_idx]
    unique_frames, start_idx, counts = np.unique(sorted_frames, return_index=True, return_counts=True)
    frame_to_indices = {
        frame: sort_idx[start:start + count]
        for frame, start, count in zip(unique_frames, start_idx, counts)
    }
    return unique_frames, frame_to_indices


def segment_by_num_locs_per_window(loc_frames: np.ndarray, min_n_locs_per_window: int,
                                   max_locs_per_segment = None,
                                   return_param_dict: bool = False) -> SegmentationResult:
    """
    Segments by collecting a minimum number of localizations per window.
    Once the threshold is met and enough locs remain, a new segment is created.
    This method ensures that each segment has at least `min_n_locs_per_window` localizations,
    while also trying to avoid creating segments that are too small at the end of the dataset.
    If `max_locs_per_segment` is set, a random subset of that size is chosen from each segment.
    This method is particularly useful for datasets with varying localization densities over time.
    Parameters:
    loc_frames (np.ndarray): Array of frame numbers for each localization.
    min_n_locs_per_window (int): Minimum number of localizations per segment.
    max_locs_per_segment (Optional[int]): Maximum number of localizations per segment. If None, all locs are used.
    return_param_dict (bool): Whether to return a dictionary of segmentation parameters.
    Returns:
    SegmentationResult: A dataclass containing segmentation results and parameters.
    """
    loc_frames = np.asarray(loc_frames, dtype=int)
    n_locs = len(loc_frames)

    if max_locs_per_segment is not None and max_locs_per_segment < 1:  # downsampling in percentage
        max_locs_per_segment = int(min_n_locs_per_window * max_locs_per_segment)

    unique_frames, frame_to_indices = _group_by_frame(loc_frames)
    loc_segments = np.full(n_locs, -1, dtype=int)  # Default to -1 for safety
    segment_counter = 0
    n_locs_in_current_segment = 0
    current_segment_indices = []
    start_frames, end_frames, locs_per_segment = [], [], []

    for i, frame in enumerate(unique_frames):
        indices = frame_to_indices[frame]
        n_locs_this_frame = len(indices)
        remaining_locs = n_locs - (len(current_segment_indices) + n_locs_this_frame + np.sum(locs_per_segment))

        # Add frame to current segment if:
        # - It fills the current segment to threshold
        # - AND there are enough locs left for another segment (or it's the last frame)
        if (n_locs_in_current_segment + n_locs_this_frame >= min_n_locs_per_window) and \
                (remaining_locs >= min_n_locs_per_window or i == len(unique_frames) - 1):
            current_segment_indices.extend(indices)
            n_locs_in_current_segment += n_locs_this_frame

            loc_segments[current_segment_indices] = segment_counter
            start_frames.append(loc_frames[current_segment_indices[0]])
            end_frames.append(loc_frames[current_segment_indices[-1]])
            locs_per_segment.append(len(current_segment_indices))

            segment_counter += 1
            current_segment_indices = []
            n_locs_in_current_segment = 0
        else:
            # Defer frame to current segment
            current_segment_indices.extend(indices)
            n_locs_in_current_segment += n_locs_this_frame

    n_segments = segment_counter
    center_frames = np.zeros(n_segments)
    loc_valid = np.zeros(n_locs, dtype=bool)

    for i in range(n_segments):
        segment_indices = np.where(loc_segments == i)[0]
        if max_locs_per_segment and len(segment_indices) > max_locs_per_segment:
            selected = np.random.choice(segment_indices, max_locs_per_segment, replace=False)
        else:
            selected = segment_indices
        loc_valid[selected] = True
        locs_per_segment[i] = len(selected)
        center_frames[i] = np.mean(loc_frames[selected])

    out_dict = None
    if return_param_dict:
        n_locs_valid = loc_valid.sum()
        out_dict = {
            "n_segments": n_segments,
            "min_n_locs_per_window": min_n_locs_per_window,
            "frames_per_window": -1,
            "start_frames": np.array(start_frames),
            "end_frames": np.array(end_frames),
            "locs_per_segment": np.array(locs_per_segment),
            "n_locs": n_locs,
            "n_locs_valid": n_locs_valid,
            "n_locs_invalid": n_locs - n_locs_valid,
            "center_frames": center_frames
        }

    return SegmentationResult(loc_segments, loc_valid, center_frames, n_segments, out_dict)


def segment_by_frame_windows(loc_frames: np.ndarray, n_frames_per_window: int,
                             max_locs_per_segment = None,
                             return_param_dict: bool = False) -> SegmentationResult:
    """
    Splits localization data into fixed-size windows of N frames.
    All localizations in those frames are grouped into one segment.
    """
    loc_frames = np.asarray(loc_frames, dtype=int)
    frames, frame_to_indices = _group_by_frame(loc_frames)
    n_locs = len(loc_frames)
    n_segments = int(np.ceil(len(frames) / n_frames_per_window))

    if max_locs_per_segment is not None and max_locs_per_segment < 1:  # downsampling in percentage
        max_locs_per_segment = n_locs / n_segments * max_locs_per_segment

    loc_segments = np.zeros(n_locs, dtype=int)
    center_frames = np.zeros(n_segments)
    loc_valid = np.ones(n_locs, dtype=bool)
    start_frames, end_frames, locs_per_segment = [], [], []

    for i in range(n_segments):
        frame_window = frames[i * n_frames_per_window:(i + 1) * n_frames_per_window]
        indices = np.concatenate([frame_to_indices[frame] for frame in frame_window if frame in frame_to_indices])
        if len(indices) == 0:
            continue
        loc_segments[indices] = i
        start_frames.append(frame_window[0])
        end_frames.append(frame_window[-1])
        center_frames[i] = np.mean(loc_frames[indices])
        locs_per_segment.append(len(indices))
        if max_locs_per_segment and len(indices) > max_locs_per_segment:
            mask = np.ones(len(indices), dtype=bool)
            mask[np.random.choice(len(indices), len(indices) - max_locs_per_segment, replace=False)] = False
            loc_valid[indices[~mask]] = False

    out_dict = None
    if return_param_dict:
        n_locs_valid = loc_valid.sum()
        out_dict = {
            "n_segments": n_segments,
            "min_n_locs_per_window": -1,
            "frames_per_window": n_frames_per_window,
            "start_frames": np.array(start_frames),
            "end_frames": np.array(end_frames),
            "locs_per_segment": np.array(locs_per_segment),
            "n_locs": n_locs,
            "n_locs_valid": n_locs_valid,
            "n_locs_invalid": n_locs - n_locs_valid,
            "center_frames": center_frames
        }

    return SegmentationResult(loc_segments, loc_valid, center_frames, n_segments, out_dict)


def segment_by_num_windows(loc_frames: np.ndarray, n_windows: int, max_locs_per_segment = None,
                           return_param_dict: bool = False) -> SegmentationResult:
    """
        Converts number of windows into an equivalent minimum locs per window,
        then calls `segment_by_num_locs_per_window`.
    """
    n_locs = len(loc_frames)
    n_locs_per_window = int(np.ceil(n_locs / n_windows))
    if max_locs_per_segment is not None and max_locs_per_segment < 1: # downsampling in percentage
        max_locs_per_segment = int(n_locs_per_window * max_locs_per_segment)
    return segment_by_num_locs_per_window(loc_frames, n_locs_per_window, max_locs_per_segment, return_param_dict)


def segmentation_wrapper(loc_frames: np.ndarray, segmentation_var: int, segmentation_mode: int = 2,
                         max_locs_per_segment = None,
                         return_param_dict: bool = False) -> SegmentationResult:
    """
        Dispatch function that selects segmentation method:
        0 → fixed number of windows
        1 → fixed number of localizations per window
        2 → fixed number of frames per window (default)
    """
    if segmentation_mode == 0:
        return segment_by_num_windows(loc_frames, segmentation_var, max_locs_per_segment, return_param_dict)
    elif segmentation_mode == 1:
        return segment_by_num_locs_per_window(loc_frames, segmentation_var, max_locs_per_segment, return_param_dict)
    else:
        return segment_by_frame_windows(loc_frames, segmentation_var, max_locs_per_segment, return_param_dict)
