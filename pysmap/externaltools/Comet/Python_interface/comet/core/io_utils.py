import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import csv

from comet.core._dialogs import ask_open_filename, ask_save_filename


def load_thunderstorm_csv(filename=None, return_essentials=True):
    """
    Load a ThunderSTORM CSV file and return the localization data.
    Parameters:
    - filename: str, path to the CSV file. If None, a file dialog will open.
    - return_essentials: bool, if True, return only essential columns (x, y, z, frame).
    Returns:
    - np.ndarray or pd.DataFrame: localization data.
    """
    if filename is None:
        filename = ask_open_filename(title="Select ThunderSTORM CSV file")

    try:
        df = pd.read_csv(filename)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        raise FileNotFoundError("Failed to read ThunderSTORM CSV file.")

    if not return_essentials:
        return df

    keys = df.columns.to_numpy()
    print(f"Opened ThunderSTORM file with keys: {keys}")

    # Initialise locs information with following convention:
    # Column 0: x coordinates
    # Column 1: y coordinates
    # Column 0: z coordinates
    # Column 0: time frames
    locs = np.zeros((len(df), 4))

    locs[:, 0] = df["x [nm]"]
    locs[:, 1] = df["y [nm]"]
    # Fill the z coordinates in 3rd column if available, leave full of 0's otherwise
    if "z [nm]" in keys:
        locs[:, 2] = df["z [nm]"]
    locs[:, 3] = df["frame"]

    return locs


# Molecule-set files store positions in pixel units and the pixel size alongside
# them. The key naming changed between format versions: older files carry a
# single general pixel size, newer ones carry separate xy and z values. The
# axis-specific keys are checked first so that a file carrying both does not get
# its z scaled by the general value.
_XY_PIXEL_SIZE_KEYS = ('xy_pixel_size_um', 'pixel_size_um', 'pixelsize_um')
_Z_PIXEL_SIZE_KEYS = ('z_pixel_size_um',)


def _first_scalar(group, keys):
    """Return the first key present in `group` as a float, or None."""
    for key in keys:
        if key in group:
            value = np.asarray(group[key]).reshape(-1)
            if value.size:
                return float(value[0])
    return None


def _read_pixel_sizes(group):
    """Read (xy, z) pixel sizes in nm from a molecule_set_data group.

    Falls back to the xy pixel size for z when the file predates the separate
    z key, which is what a single general pixel size means.
    """
    xy_um = _first_scalar(group, _XY_PIXEL_SIZE_KEYS)
    if xy_um is None:
        raise KeyError(
            "No pixel size found in molecule_set_data; expected one of {}.".format(
                ", ".join(_XY_PIXEL_SIZE_KEYS))
        )
    z_um = _first_scalar(group, _Z_PIXEL_SIZE_KEYS)
    if z_um is None:
        z_um = xy_um
    return xy_um * 1e3, z_um * 1e3


def load_normal_molecule_set(filename=None, sanity_check=False, photon_bandpass=(None, None)):
    """
    Load a normal molecule set from an HDF5 file.
    Parameters:
    - filename: str, path to the HDF5 file. If None, a file dialog will open.
    - sanity_check: bool, if True, display a scatter plot of the XY positions.
    - photon_bandpass: tuple, (min_photons, max_photons) to filter localizations by photon count.
    Returns:
    - np.ndarray: localization data with columns [x_nm, y_nm, z_nm, frame].
    """
    if filename is None:
        filename = ask_open_filename(title="Select HDF5 dataset")

    f = h5py.File(filename, 'r')
    pixelsize_nm, pixelsize_z_nm = _read_pixel_sizes(f['molecule_set_data'])

    datatable = f['molecule_set_data']['datatable']
    x_pos = np.asarray(datatable['X_POS_PIXELS']) * pixelsize_nm
    y_pos = np.asarray(datatable['Y_POS_PIXELS']) * pixelsize_nm
    try:
        z_pos = np.asarray(datatable['Z_POS_PIXELS']) * pixelsize_z_nm
    except Exception:
        print("z coordinates not found in molecule set, proceeding with 2D data")
        z_pos = np.zeros_like(x_pos)
    frames = np.asarray(datatable['FRAME_NUMBER'])

    locs = np.stack([x_pos, y_pos, z_pos, frames], axis=1)

    if photon_bandpass[0] is not None and photon_bandpass[1] is not None:
        photons = np.asarray(datatable['PHOTONS'])
        mask = (photons > photon_bandpass[0]) & (photons < photon_bandpass[1])
        locs = locs[mask]

    if sanity_check:
        plt.figure()
        plt.scatter(locs[:, 0], locs[:, 1], alpha=0.01)
        plt.title("Sanity Check: XY Projection")
        plt.show()

    return locs


def load_simulation_dataset_and_gt_drift(filename=None, display=False, remove_loc_prec=False):
    if filename is None:
        filename = ask_open_filename(title="Select simulation dataset")
    f = h5py.File(filename)
    frames = np.asarray(f['frame_number'], dtype=int)
    drift = np.asarray(f['sample_drift']['drift_data']) * 1E3
    if not remove_loc_prec:
        x = np.asarray(f['x_coords']) * 1E3
        y = np.asarray(f['y_coords']) * 1E3
        z = np.asarray(f['z_coords']) * 1E3
    else:
        label_sites_nm = np.asarray(f['label_sites']['site_data']) * 1E3
        label_site_idc = np.asarray(f['label_site_index'])
        x = label_sites_nm[0, label_site_idc]
        y = label_sites_nm[1, label_site_idc]
        z = label_sites_nm[2, label_site_idc]
        x += drift[0, frames]
        y += drift[1, frames]
        z += drift[2, frames]

    f.close()
    print(f"Imported {len(x)} coords on {frames.max() + 1} frames with a max drift of {np.max(np.abs(drift))} nm")
    if display:
        plt.figure()
        plt.scatter(x[::10], y[::10], alpha=0.01)
        plt.figure()
        plt.plot(drift[0, :])
        plt.plot(drift[1, :])
        plt.plot(drift[2, :])
        plt.show()
    dataset = np.zeros((len(x), 4))
    dataset[:, 0] = x
    dataset[:, 1] = y
    dataset[:, 2] = z
    dataset[:, 3] = frames
    drift = np.swapaxes(drift, 0, 1)
    return dataset, drift


def save_drift_correction_details(savename, drift_est, drift_est_intp, frames_intp,
                                  segmentation_result, calc_time, initial_sigma_nm, target_sigma_nm, gt_drift=None):
    """
    Save drift correction details to an HDF5 file.
    Parameters:
    - savename: str, path to save the HDF5 file.
    - drift_est: np.ndarray, estimated drift per segment in nm.
    - drift_est_intp: np.ndarray, interpolated drift estimates in nm.
    - frames_intp: np.ndarray, frames corresponding to the interpolated drift estimates.
    - segmentation_result: SegmentationResult, result of the segmentation process.
    - calc_time: float, time taken for the drift calculation in seconds.
    - initial_sigma_nm: float, initial sigma used in drift estimation.
    - target_sigma_nm: float, target sigma used in drift estimation.
    - gt_drift: np.ndarray or None, ground truth drift if available.
    """
    f = h5py.File(savename, 'a')

    f.create_group("drift")
    f['drift']['center_frames'] = segmentation_result.center_frames
    f['drift']['drift_per_segment_nm'] = drift_est
    f['drift']['drift_interpolated_nm'] = drift_est_intp
    f['drift']['frames_interpolated'] = frames_intp
    if gt_drift is not None:
        f['drift']['gt_drift'] = gt_drift
    f.create_group('parameters')
    f['parameters']['calculation_time_s'] = calc_time
    for key, value in segmentation_result.out_dict.items():
        f['parameters'][key] = value
    f['parameters']['initial_sigma_nm'] = initial_sigma_nm
    f['parameters']['target_sigma_nm'] = target_sigma_nm
    f.close()


def save_dataset_as_ms_h5(storm_coordinates, frames, pixelsize_nm, pixelsize_z_nm=None, amplitudes=None,
                                 uncertainty_x=None, uncertainty_y=None, uncertainty_z=None, frame_shape=None,
                                 filename=None, extra_dict=None):
    """
    Save dataset in Molecule Set format used by Daxview.
    Parameters:
    - storm_coordinates: np.ndarray, shape (N, 2) or (N, 3), coordinates in nm.
    - frames: np.ndarray, shape (N,), frame numbers.
    - pixelsize_nm: float, pixel size in nm.
    - pixelsize_z_nm: float or None, pixel size in z in nm. If None, set equal to pixelsize_nm.
    - amplitudes: np.ndarray or None, shape (N,), photon counts. If None, set to 1000.
    - uncertainty_x: np.ndarray or None, shape (N,), uncertainty in x in nm. If None, set to 10 nm.
    - uncertainty_y: np.ndarray or None, shape (N,), uncertainty in y in nm. If None, set to 10 nm.
    - uncertainty_z: np.ndarray or None, shape (N,), uncertainty in z in nm. If None, set to 10 nm.
    - frame_shape: tuple or None, shape of the frames (height, width). If None, set to (128, 128).
    - filename: str or None, path to save the HDF5 file. If None, a file dialog will open.
    - extra_dict: dict or None, additional key-value pairs to save in the file.
    """
    if filename is None:
        filename = ask_save_filename(defaultextension=".h5", filetypes=[("hdf5 files", "*.h5")])
    f = h5py.File(filename, 'w')

    headers = ["X_POS_PIXELS", "Y_POS_PIXELS", "Z_POS_PIXELS", "PRECISION_XY_PIXELS", "PRECISION_Z_PIXELS",
               "PHOTONS", "CHANNEL", "FRAME_NUMBER", "INDEX"]
    dtypes = [np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.int32, np.int32, np.int32]
    compound_dtype = np.dtype([(headers[i], dtypes[i]) for i in range(len(headers))])

    if pixelsize_z_nm is None:
        pixelsize_z_nm = pixelsize_nm

    structured_array = np.zeros(len(frames), dtype=compound_dtype)
    # A molecule set stores positions in pixel units alongside the pixel size,
    # so every coordinate is divided by the pixel size of its own axis. The
    # first column is X; there is no axis swap.
    structured_array['X_POS_PIXELS'] = storm_coordinates[:, 0] / pixelsize_nm
    structured_array['Y_POS_PIXELS'] = storm_coordinates[:, 1] / pixelsize_nm
    if storm_coordinates.shape[1] > 2:
        structured_array['Z_POS_PIXELS'] = storm_coordinates[:, 2] / pixelsize_z_nm
    else:
        structured_array['Z_POS_PIXELS'] = np.zeros(len(frames))
    if amplitudes is not None:
        structured_array['PHOTONS'] = amplitudes
    else:
        structured_array['PHOTONS'] = np.ones(len(frames)) * 1E3
    if uncertainty_x is not None:
        structured_array['PRECISION_XY_PIXELS'] = np.sqrt(uncertainty_x ** 2 + uncertainty_y ** 2) / pixelsize_nm
    else:
        structured_array['PRECISION_XY_PIXELS'] = np.ones(len(frames)) / 10
    if uncertainty_z is not None:
        structured_array['PRECISION_Z_PIXELS'] = uncertainty_z / pixelsize_z_nm
    else:
        structured_array['PRECISION_Z_PIXELS'] = np.ones(len(frames)) / 10
    structured_array['CHANNEL'] = np.zeros(len(frames))
    structured_array['FRAME_NUMBER'] = frames
    structured_array['INDEX'] = np.arange(len(frames))

    f['daxview_file_type'] = np.array('DAXVIEW_STORM_DATA', dtype=h5py.string_dtype('ascii', 19))
    f['object_class_name'] = np.array('DV_STORM_4PI_MOLECULE_SET', dtype=h5py.string_dtype('ascii', 26))
    f.create_dataset('save_format_version', data=np.ones(1, dtype=np.float32) * 1.2)
    f.create_group('molecule_set_data')
    f['molecule_set_data'].create_dataset('datatable', data=structured_array)
    f['molecule_set_data']['flag_dcorr_applied'] = 1
    f['molecule_set_data']['flag_tform_applied'] = 0
    f['molecule_set_data']['n_keys'] = 6
    f['molecule_set_data']['keynames'] = np.asarray(["z_pixel_size_um", "xy_pixel_size_um", "version", "n_molecules",
                                                     "flag_dcorr_applied", "xpixels", "stat_t_matrix", "ypixels",
                                                     "flag_tform_applied", "datatable"],
                                                    dtype=h5py.string_dtype('ascii', 19))
    f['molecule_set_data']['n_molecules'] = len(frames)
    f['molecule_set_data']['stat_t_matrix'] = np.diag([1, 1, 1, 1])
    f['molecule_set_data']['version'] = 1
    if frame_shape is None:
        f['molecule_set_data']['xpixels'] = 128
        f['molecule_set_data']['ypixels'] = 128
    else:
        f['molecule_set_data']['xpixels'] = frame_shape[0]
        f['molecule_set_data']['ypixels'] = frame_shape[1]
    f['molecule_set_data']['xy_pixel_size_um'] = pixelsize_nm / 1E3
    f['molecule_set_data']['z_pixel_size_um'] = pixelsize_z_nm / 1E3

    if extra_dict is not None:
        group = f.create_group('extra_data')
        for key, value in extra_dict.items():
            group[key] = value
    f.close()


def correct_and_save_thunderstorm_csv(drift_interp_with_frames_nm, filename=None, savename=None):
    """
    Load a ThunderSTORM CSV file, apply drift correction, and save the corrected dataset.
    Parameters:
    - drift_interp_with_frames_nm: np.ndarray, shape (num_frames, 3) or (num_frames, 4)
      for 2D or 3D datasets + frames respectively, drift values in nm.
    - filename: str or None, path to the input CSV file. If None, a file dialog will open.
    - savename: str or None, path to save the corrected CSV file. If None, a file dialog will open.

    """
    if filename is None:
        filename = ask_open_filename(title="Select ThunderSTORM CSV file to correct")

    df = load_thunderstorm_csv(filename, return_essentials=False)
    frames = df["frame"].to_numpy().astype(int)
    n_drift_frames = drift_interp_with_frames_nm.shape[0]
    if frames.max() >= n_drift_frames:
        raise ValueError(
            f"CSV contains frame index {frames.max()} but the drift table only covers "
            f"frames 0..{n_drift_frames - 1}. The drift was probably estimated from a "
            f"different dataset than the one being corrected."
        )
    df['x [nm]'] -= drift_interp_with_frames_nm[frames, 0]
    df['y [nm]'] -= drift_interp_with_frames_nm[frames, 1]
    if "z [nm]" in df.columns:
        df['z [nm]'] -= drift_interp_with_frames_nm[frames, 2]
    if savename is None:
        savename = ask_save_filename(title="Save corrected ThunderSTORM CSV", defaultextension=".csv")
    df.to_csv(savename, index=False, quoting=csv.QUOTE_NONE, float_format='%.3f')


def save_dataset_as_thunderstorm_csv(dataset, savename=None):
    """
    Save dataset in ThunderSTORM CSV format.
    Parameters:
    - dataset: np.ndarray, shape (N, 4), localization data with columns [x_nm, y_nm, z_nm, frame].
    - savename: str or None, path to save the CSV file. If None, a file dialog will open.
    """
    if savename is None:
        savename = ask_save_filename(title="Save as ThunderSTORM CSV", defaultextension=".csv")

    df = pd.DataFrame({
        "frame": dataset[:, 3].astype(int),
        "x [nm]": dataset[:, 0],
        "y [nm]": dataset[:, 1],
    })

    header = ['"frame"', '"x [nm]"', '"y [nm]"']
    # Add z coordinates if present
    if (np.any(dataset[:, 2] != np.zeros(len(dataset)))):
        df["z [nm]"] = dataset[:, 2]
        header.append('"z [nm]"')

    df.to_csv(savename, index=True,
              header=header, index_label='"id"',
              mode='w', quoting=csv.QUOTE_NONE, 
              float_format = '%.3f') # precision of 3 


def save_dataset_custom(dataset, savename=None, format_hint=None):
    """
    Placeholder function. Extend this to support custom save formats.
    Use 'format_hint' to route specific formats.
    Parameters:
    - dataset: np.ndarray, localization data.
    - savename: str or None, path to save the file. If None, a file dialog will open.
    - format_hint: str or None, hint for the desired format.
    """
    raise NotImplementedError("Custom save format handler not implemented. Use HDF5 or ThunderSTORM CSV.")
