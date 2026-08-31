import h5py
import numpy as np
import glob
import matplotlib.pyplot as plt
from comet.core.drift_optimizer import comet_run_kd
from comet.core._dialogs import ask_directory

from comet.core.io_utils import save_dataset_as_thunderstorm_csv
from comet.utilities.post_analysis_utilities import analyse_folder_of_h5_drift_summary_files

rendering_pixelsize_nm = 10

def analyse_ai_minflux_npy_files_for_drift(folder=None, max_drift_nm=100, traces_per_time_window=50, save_result_as_csv=False, sanity_check=False):
    """
        Analyze .npy files from AI Minflux for drift correction using COMET. --> Old .npy format (b4 2025)
        Parameters:
        - folder: Directory containing .npy files. If None, a dialog will prompt for
            folder selection.
        - max_drift_nm: Maximum expected drift in nanometers.
        - traces_per_time_window: Number of traces per time window for segmentation.
        - save_result_as_csv: If True, saves the drift-corrected dataset as a
            ThunderSTORM-compatible CSV file.
        - sanity_check: If True, displays intermediate histograms for verification.
    """
    if folder is None:
        folder = ask_directory(title="Select folder with .npy files")
    files = glob.glob(folder + "/*.npy")
    for file in files:
        data = np.load(file)
        locs = data['itr']['loc'][:, -1] * 1E9  # convert to nm
        tid = data['tid']
        utid, utid_idc, utid_counts = np.unique(tid, return_index=True, return_counts=True)
        print(f"Processing {len(locs)} localizations in {len(utid)} traces from file: {file}.")
        mean_loc_per_trace =[]
        frame_to_trace_id = []
        for i, id in enumerate(utid):
            #if utid_counts[i] > n_locs_per_trace:
            idc_start = utid_idc[i]
            idc_end = idc_start + utid_counts[i]
            locs_id = locs[idc_start:idc_end]
            mean_loc_per_trace.append(np.mean(locs_id, axis=0))
            frame_to_trace_id.append(id)

        dataset = np.zeros((len(mean_loc_per_trace), 4), dtype=np.float32)
        dataset[:, :3] = mean_loc_per_trace  # x,y,z
        dataset[:, 3] = np.arange(len(mean_loc_per_trace))    # frame

        for i in range(3): # hsift to origin for convenience
            dataset[:, i] -= np.min(dataset[:, i])
        if sanity_check:
            plt.figure()
            plt.hist2d(dataset[:, 0], dataset[:, 1],
                       bins=[np.arange(np.max(dataset[:, 0])//rendering_pixelsize_nm)*rendering_pixelsize_nm,
                       np.arange(np.max(dataset[:, 1])//rendering_pixelsize_nm)*rendering_pixelsize_nm])
            plt.show()

        drift_nm, dataset_corr = comet_run_kd(
            dataset=dataset.copy(),
            segmentation_mode=2,
            segmentation_var=traces_per_time_window,
            initial_sigma_nm=max_drift_nm//3,
            max_drift_nm=max_drift_nm,
            target_sigma_nm=10,
            drift_max_bound_factor=2,
            boxcar_width=3,
            return_corrected_locs=True,
            save_correction_details=True,
            interpolation_method='cubic',
            mode="torch",
            display=True
        )

        if sanity_check:
            fig, ax = plt.subplots(1,2, sharex=True, sharey=True)
            ax[0].hist2d(dataset[:, 0], dataset[:, 1],
                       bins=[np.arange(np.max(dataset[:, 0])//rendering_pixelsize_nm)*rendering_pixelsize_nm,
                       np.arange(np.max(dataset[:, 1])//rendering_pixelsize_nm)*rendering_pixelsize_nm])
            for i in range(3):
                dataset_corr[:, i] -= np.min(dataset_corr[:, i])
            ax[1].hist2d(dataset_corr[:, 0], dataset_corr[:, 1],
                         bins=[np.arange(np.max(dataset_corr[:, 0]) // rendering_pixelsize_nm) * rendering_pixelsize_nm,
                               np.arange(np.max(dataset_corr[:, 1]) // rendering_pixelsize_nm) * rendering_pixelsize_nm])
            plt.show()


        locs_corr = locs.copy()

        for i in range(len(utid)):
            idc = np.arange(utid_counts[i])+utid_idc[i]
            locs_corr[idc] -= drift_nm[i, :-1]

        if save_result_as_csv:
            save_dataset_as_thunderstorm_csv(dataset)

        for i in range(3):
            locs[:, i] -= np.min(locs[:, i])
            locs_corr[:, i] -= np.min(locs_corr[:, i])

        fig, ax = plt.subplots(1,2, sharex=True, sharey=True)
        ax[0].hist2d(locs[:, 0], locs[:, 1],
                   bins=[np.arange(np.max(locs)//rendering_pixelsize_nm)*rendering_pixelsize_nm,
                   np.arange(np.max(locs[:, 1])//rendering_pixelsize_nm)*rendering_pixelsize_nm])

        ax[1].hist2d(locs_corr[:, 0], locs_corr[:, 1],
                     bins=[np.arange(np.max(locs_corr[:, 0]) // rendering_pixelsize_nm) * rendering_pixelsize_nm,
                           np.arange(np.max(locs_corr[:, 1]) // rendering_pixelsize_nm) * rendering_pixelsize_nm])
        plt.show()







if __name__ == "__main__":
    #folder = ask_directory(title="Select folder with .npy files")
    #analyse_ai_minflux_npy_files_for_drift(folder, max_drift_nm=300, traces_per_time_window=50, sanity_check=False, save_result_as_csv=True)

    analyse_folder_of_h5_drift_summary_files()

