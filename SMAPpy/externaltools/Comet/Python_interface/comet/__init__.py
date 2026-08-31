"""COMET - Cost-function Optimized Maximal Overlap Drift EsTimation.

Drift correction for single-molecule localization microscopy (SMLM) datasets.

Typical use::

    import numpy as np
    from comet import comet_run_kd, load_thunderstorm_csv

    locs = load_thunderstorm_csv("my_data.csv")       # (N, 4): x, y, z, frame
    drift, corrected = comet_run_kd(
        dataset=locs,
        segmentation_mode=2,      # frames per time window
        segmentation_var=50,
        max_drift_nm=100,
        return_corrected_locs=True,
    )

The compute backend defaults to ``"cuda"``; pass ``mode="cpu"`` or
``mode="torch"`` to override, or use :func:`comet.best_backend` to pick
whatever is fastest on the current machine.
"""

from comet._version import __version__

from comet.core.backends import (
    available_backends,
    best_backend,
    cuda_available,
    describe_backends,
    torch_available,
)
from comet.core.drift_optimizer import comet_run_kd
from comet.core.io_utils import (
    correct_and_save_thunderstorm_csv,
    load_normal_molecule_set,
    load_thunderstorm_csv,
    save_dataset_as_ms_h5,
    save_dataset_as_thunderstorm_csv,
)
from comet.core.segmenter import segmentation_wrapper

__all__ = [
    "__version__",
    # pipeline
    "comet_run_kd",
    "segmentation_wrapper",
    # io
    "load_thunderstorm_csv",
    "load_normal_molecule_set",
    "save_dataset_as_thunderstorm_csv",
    "save_dataset_as_ms_h5",
    "correct_and_save_thunderstorm_csv",
    # backend introspection
    "available_backends",
    "best_backend",
    "cuda_available",
    "torch_available",
    "describe_backends",
]
