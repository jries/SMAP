"""Open an acquisition, report what the metadata gives us, read a few frames."""
import sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from smappy.io.tiff import open_stack, metadata_from_stack
from smappy.io.cameras_mat import CameraPresets

src = open_stack(sys.argv[1])
presets = CameraPresets.load(sys.argv[2]) if len(sys.argv) > 2 else None
print(f"files      : {[f.name for f in src.files]}")
print(f"frames     : {src.n_frames} written, {src.n_frames_declared} declared by MM")
print(f"frame      : {src.shape} {src.dtype}")
if presets: print(f"presets    : {presets.describe(src.mm_metadata)}")
cam = metadata_from_stack(src, presets)
print(f"metadata   : {cam}")
print(f"             exposure={cam.exposure_ms} ms  excess_noise={cam.excess_noise}")
try:
    print(f"             adu->photons factor = {cam.adu_to_photons:.5f}")
except Exception as e:
    print("             ", e)
try:
    cam.require()
except ValueError as e:
    print(f"MISSING    : {e}")

t = time.time()
n = 0
for start, block in src.frames(chunk=50, stop=200):
    n += len(block)
    if start == 0:
        print(f"\nfirst block: start={start} shape={block.shape} "
              f"min={block.min()} max={block.max()} median={np.median(block):.0f}")
print(f"read {n} frames in {time.time()-t:.2f} s")
last = src.frame(src.n_frames - 1)
print(f"last frame : median={np.median(last):.0f} (non-zero => not padding)")

# optional third argument: a YAML file with user overrides
if len(sys.argv) > 3:
    from smappy.io.tiff import camera_metadata
    cam2 = camera_metadata(src, presets, sys.argv[3])
    print(f"\nwith config : {cam2}")
    print(f"              complete, adu->photons = {cam2.adu_to_photons:.5f}")
