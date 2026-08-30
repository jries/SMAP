import numpy as np
from scipy.spatial import cKDTree
import h5py


def pair_indices_kdtree(coordinates, distance):
    """
    Find all pairs of points within a certain distance using a KD-tree.
    Parameters:
    - coordinates: np.ndarray of shape (N, D) where N is the number of points and D is the dimensionality.
    - distance: float, the maximum distance to consider points as a pair.
    Returns:
    - idx1: np.ndarray of shape (M,), indices of the first point in each pair.
    - idx2: np.ndarray of shape (M,), indices of the second point in each pair.
    """
    tree = cKDTree(coordinates)
    try:
        pairs = tree.query_pairs(r=distance, output_type='ndarray')
    except MemoryError:
        print("[pair_indices_kdtree] MemoryError encountered")
        # Typed empties so callers can use .size/.astype without special-casing.
        empty = np.empty(0, dtype=np.int32)
        return empty, empty, False
    return np.ascontiguousarray(pairs[:, 0], dtype=np.int32), np.ascontiguousarray(pairs[:, 1], dtype=np.int32), True

def pair_indices_kdtree_legacy(coordinates, distance):
    """
    Find all pairs of points within a certain distance using a KD-tree.
    Parameters:
    - coordinates: np.ndarray of shape (N, D) where N is the number of points and D is the dimensionality.
    - distance: float, the maximum distance to consider points as a pair.
    Returns:
    - idx1: np.ndarray of shape (M,), indices of the first point in each pair.
    - idx2: np.ndarray of shape (M,), indices of the second point in each pair.
    """
    tree = cKDTree(coordinates)
    while True:
        try:
            pairs = tree.query_pairs(r=distance, output_type='ndarray')
            break
        except MemoryError:
            distance *= 0.8
            print(f"[pair_indices_kdtree] Reducing distance to {distance:.2f} due to memory error.")

    print(f"[pair_indices_kdtree] Found {len(pairs):,} pairs")
    return np.ascontiguousarray(pairs[:, 0], dtype=np.int32), np.ascontiguousarray(pairs[:, 1], dtype=np.int32)


def pair_indices_kdtree_full_to_file_recursion(coordinates, distance, filename=None, split_dimension=0, indices=None):
    """
    Recursively find all pairs of points within a certain distance using a KD-tree,
    and save the results to an HDF5 file to avoid memory issues.
    Parameters:
    - coordinates: np.ndarray of shape (N, D) where N is the number of points and D is the dimensionality.
    - distance: float, the maximum distance to consider points as a pair.
    - filename: str, path to the HDF5 file where results will be saved.
    - split_dimension: int, the dimension along which to split the data when a MemoryError occurs.
    - indices: np.ndarray of shape (N,), original indices of the points in the full dataset.
    Returns:
    - None (results are saved to the specified HDF5 file)
    """
    if indices is None:
        indices = np.arange(len(coordinates))

    try:
        tree = cKDTree(coordinates)
        pairs = tree.query_pairs(r=distance, output_type='ndarray')
        global_pairs = np.stack([indices[pairs[:, 0]], indices[pairs[:, 1]]], axis=1)

        with h5py.File(filename, 'a') as f:
            grp_name = f"pair_indices_{len(f.keys()):04d}"
            f.create_dataset(grp_name, data=global_pairs, compression="gzip")
            print(f"[pair_indices_kdtree_full_to_file_recursion] Saved {len(global_pairs):,} pairs to '{grp_name}'")

    except MemoryError:
        print("[pair_indices_kdtree_full_to_file_recursion] MemoryError - splitting...")
        median_val = np.median(coordinates[:, split_dimension])
        left_mask = coordinates[:, split_dimension] < median_val
        right_mask = ~left_mask
        dim_next = (split_dimension + 1) % coordinates.shape[1]

        pair_indices_kdtree_full_to_file_recursion(
            coordinates[left_mask], distance, filename,
            split_dimension=dim_next, indices=indices[left_mask]
        )
        pair_indices_kdtree_full_to_file_recursion(
            coordinates[right_mask], distance, filename,
            split_dimension=dim_next, indices=indices[right_mask]
        )


def pair_indices_lex_floor_asymmetric(coordinates, distance):
    coordinates = coordinates.copy()  # the loop below shifts coordinates in place
    for i in range(len(coordinates[0])):
        coordinates[:, i] -= np.min(coordinates[:, i])
    coordinates = np.array(np.floor(coordinates / distance), dtype=int)
    coordinates = np.array(list(map(tuple, coordinates)))

    sort_indices = np.lexsort(coordinates.T)

    if coordinates.shape[1] == 2:
        tmp = coordinates[:, 0].copy()
        coordinates[:, 0] = coordinates[:, 1]
        coordinates[:, 1] = tmp
    else:
        tmp = coordinates[:, 0].copy()
        coordinates[:, 0] = coordinates[:, 1]
        coordinates[:, 1] = tmp
        tmp = coordinates[:, 0].copy()
        coordinates[:, 0] = coordinates[:, 2]
        coordinates[:, 2] = tmp
        tmp = coordinates[:, 2].copy()
        coordinates[:, 2] = coordinates[:, 1]
        coordinates[:, 1] = tmp
    # get the unique tuples and their counts
    unique_tuples, counts = np.unique(coordinates[sort_indices], axis=0, return_counts=True)
    # get the indices of the similar tuples
    similar_indices = np.split(sort_indices, np.cumsum(counts[:-1]))
    idx_i = []
    idx_j = []
    pair_idc_estimate = 0
    for i in range(len(similar_indices)):
        n_elements = len(similar_indices[i])
        pair_idc_estimate += 0.5 * n_elements * (n_elements + 1)
    if pair_idc_estimate > 5E8:  # default 5E8
        print(pair_idc_estimate)
        print("memory Error")
        #raise MemoryError(f"To many Pair indices for RAM ({pair_idc_estimate} pairs)")
        return [], [], False
    else:
        try:
            for i in range(len(similar_indices)):
                if len(similar_indices[i]) > 1:
                    idc = similar_indices[i]
                    for j in range(len(idc) - 1):
                        tmp_n_entries = (len(idc) - j) - 1
                        idx_i.append(np.repeat(np.int32(idc[j]), tmp_n_entries))
                        idx_j.append(idc[np.arange(tmp_n_entries, dtype=np.int32) + (j + 1)])
            return np.ascontiguousarray(np.concatenate(idx_i).ravel(), dtype=np.int32), \
                   np.ascontiguousarray(np.concatenate(idx_j).ravel(), dtype=np.int32), True
        except Exception as e:
            print(e)
            return [], [], False

def estimate_pairs(coordinates, distance):
  coordinates = coordinates.copy()  # the loop below shifts coordinates in place
  for i in range(len(coordinates[0])):
    coordinates[:, i] -= np.min(coordinates[:, i])
  coordinates = np.array(np.floor(coordinates / distance), dtype=int)
  coordinates = np.array(list(map(tuple, coordinates)))
  sort_indices = np.lexsort(coordinates.T)# get the unique tuples and their counts
  unique_tuples, counts = np.unique(coordinates[sort_indices], axis=0, return_counts=True)
  # get the indices of the similar tuples
  similar_indices = np.split(sort_indices, np.cumsum(counts[:-1]))
  idx_i = []
  idx_j = []
  pair_idc_estimate = 0
  for i in range(len(similar_indices)):
    n_elements = len(similar_indices[i])
    pair_idc_estimate += n_elements * (n_elements - 1)
  rounded = round(pair_idc_estimate,-4)
  return rounded


