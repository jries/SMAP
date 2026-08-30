import glob
import os

import numpy as np
import h5py

from comet.core._dialogs import ask_directory, ask_open_filename, ask_save_filename


def combine_two_drift_estimates_from_correction_details_files(file_1, file_2, savename=None, sanity_check=False):
    if savename is None:
        savename = ask_save_filename(defaultextension=".h5", title="Save combined drift as...")
    with h5py.File(file_1, 'r') as f1, h5py.File(file_2, 'r') as f2, h5py.File(savename, 'a') as fout:

        assert np.array_equal(f1['drift']['frames_interpolated'][:], f2['drift']['frames_interpolated'][:]), "Frame numbers do not match between the two files."
        fout.require_group('drift_correction_1')
        for key in f1:
            f1.copy(key, fout['drift_correction_1'])
        fout.require_group('drift_correction_2')
        for key in f2:
            f2.copy(key, fout['drift_correction_2'])

        drift_1 = f1['drift']['drift_interpolated_nm'][:]
        drift_2 = f2['drift']['drift_interpolated_nm'][:]

        combined_drift = (drift_1 + drift_2)
        fout.require_group('drift')
        fout['drift'].create_dataset('frames_interpolated', data=f1['drift']['frames_interpolated'][:])
        fout['drift'].create_dataset('drift_interpolated_nm', data=combined_drift)

        combined_drift_with_frames = np.zeros((combined_drift.shape[0], combined_drift.shape[1] + 1), dtype=combined_drift.dtype)
        combined_drift_with_frames[:, -1] = f1['drift']['frames_interpolated'][:]
        combined_drift_with_frames[:, :-1] = combined_drift

        if sanity_check:
            import matplotlib.pyplot as plt
            frames = combined_drift_with_frames[:, -1]
            plt.figure()
            plt.plot(frames, combined_drift_with_frames[:, 0], label='X Drift')
            plt.plot(frames, combined_drift_with_frames[:, 1], label='Y Drift')
            plt.plot(frames, combined_drift_with_frames[:, 2], label='Z Drift')
            plt.xlabel('Frame')
            plt.ylabel('Drift (nm)')
            plt.title('Combined Drift Estimate')
            plt.legend()
            plt.show()

    return combined_drift_with_frames


def analyse_folder_of_h5_drift_summary_files(folder=None):
    """Overlay the X/Y drift traces of every .h5 drift summary file in a folder."""
    import matplotlib.pyplot as plt

    if folder is None:
        folder = ask_directory(title="Select folder with .h5 drift summary files")
    files = sorted(glob.glob(os.path.join(folder, "*.h5")))

    plt.figure()
    for file in files:
        label = os.path.basename(file)
        with h5py.File(file, 'r') as f:
            drift = f['drift']['drift_interpolated_nm'][:]

            plt.plot(drift[:, 0] - np.mean(drift[:, 0]), label=f"{label} X")
            plt.plot(drift[:, 1] - np.mean(drift[:, 1]), label=f"{label} Y")
    plt.xlabel("Frame")
    plt.ylabel("Drift (nm)")
    plt.title("Drift Estimates Comparison")
    plt.legend()
    plt.show()


def create_nice_plot_from_x_and_y_drift(file=None, min_frame=0, zoom=False):
    """Publication-style X/Y drift figure from a single correction details file."""
    import matplotlib
    import matplotlib.pyplot as plt

    if file is None:
        file = ask_open_filename(title="Select correction details file...", defaultextension=".h5")
    with h5py.File(file, 'r') as f:
        drift = f['drift']['drift_interpolated_nm'][:]
        frames = f['drift']['frames_interpolated'][:]

        start_idx = np.where(frames >= min_frame)[0][0]
        drift = drift[start_idx:]
        frames = frames[start_idx:]

        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

        mm = 2 / 25.4
        matplotlib.rc('font', **{'family': 'Calibri', 'size': 10})

        plt.subplots(figsize=(43 * mm, 28 * mm))
        plt.plot(frames, drift[:, 0] - np.mean(drift[:, 0]), label='X Drift')
        plt.plot(frames, drift[:, 1] - np.mean(drift[:, 1]), label='Y Drift')
        if not zoom:
            plt.xlabel('Frame')
            plt.ylabel('Drift (nm)')
            plt.title('Drift Estimate')
            plt.legend()

        plt.show()


def compare_two_drift_estimates_from_correction_details_files(file_1, file_2):
    with h5py.File(file_1, 'r') as f1, h5py.File(file_2, 'r') as f2:
        assert np.array_equal(f1['drift']['frames_interpolated'][:], f2['drift']['frames_interpolated'][:]), "Frame numbers do not match between the two files."

        drift_1 = f1['drift']['drift_interpolated_nm'][:]
        drift_2 = f2['drift']['drift_interpolated_nm'][:]

        frames = f1['drift']['frames_interpolated'][:]

        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(frames, drift_1[:, 0], label='Drift 1 X')
        plt.plot(frames, drift_1[:, 1], label='Drift 1 Y')
        plt.plot(frames, drift_1[:, 2], label='Drift 1 Z')
        plt.plot(frames, drift_2[:, 0], '--', label='Drift 2 X')
        plt.plot(frames, drift_2[:, 1], '--', label='Drift 2 Y')
        plt.plot(frames, drift_2[:, 2], '--', label='Drift 2 Z')
        plt.xlabel('Frame')
        plt.ylabel('Drift (nm)')
        plt.title('Comparison of Two Drift Estimates')
        plt.legend()

        print(f"Diff X: {np.std(drift_1[3500:, 0] - drift_2[3500:, 0]):.2f} nm")
        print(f"Diff Y: {np.std(drift_1[3500:, 1] - drift_2[3500:, 1]):.2f} nm")
        print(f"Diff Z: {np.std(drift_1[3500:, 2] - drift_2[3500:, 2]):.2f} nm")
        plt.show()

if __name__ == "__main__":
    file_1 = ask_open_filename(title="Select first correction details file...", defaultextension=".h5")
    file_2 = ask_open_filename(title="Select second correction details file...", defaultextension=".h5")

    # combine_two_drift_estimates_from_correction_details_files(file_1, file_2, sanity_check=True)

    compare_two_drift_estimates_from_correction_details_files(file_1, file_2)