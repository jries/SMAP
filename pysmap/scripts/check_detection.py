"""Run camera conversion -> filter -> peak finding -> ROI cutting on real data.

Writes a PNG showing the photon image, the filtered image with the detected
candidates, and a montage of the first ROIs.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from smapfit.io.tiff import open_stack, camera_metadata
from smapfit.io.cameras_mat import CameraPresets
from smapfit.camera import to_photons
from smapfit.detect import DoGFilter, GaussFilter, DynamicCutoff, AbsoluteCutoff, PeakFinder
from smapfit.roi import cut_rois

data, cameras, config = sys.argv[1], sys.argv[2], sys.argv[3]
mode = sys.argv[4] if len(sys.argv) > 4 else "dog"
cutoff_value = float(sys.argv[5]) if len(sys.argv) > 5 else 1.7
nframes = int(sys.argv[6]) if len(sys.argv) > 6 else 20
roisize = int(sys.argv[7]) if len(sys.argv) > 7 else 13

src = open_stack(data)
cam = camera_metadata(src, CameraPresets.load(cameras), config)
print(f"{cam}\n")

flt = DoGFilter(1.2) if mode == "dog" else GaussFilter(1.2)
cutoff = DynamicCutoff(cutoff_value) if cutoff_value < 20 else AbsoluteCutoff(cutoff_value)
finder = PeakFinder(flt, cutoff)
print(f"{finder}, roisize={roisize}")

start, block = next(src.frames(chunk=nframes))
photons = to_photons(block, cam)
cands, filtered = finder(photons, first_frame=start)
rois = cut_rois(photons, cands, roisize, first_frame=start)

per_frame = np.bincount(cands.frame - start, minlength=nframes)
print(f"\nphotons  : median={np.median(photons):.1f} max={photons.max():.0f}")
print(f"filtered : median={np.median(filtered):.2f} p99.9={np.percentile(filtered, 99.9):.1f}")
print(f"cutoff   : {cutoff(filtered[0][filtered[0] > np.percentile(filtered[0], 99)]):.2f} (frame 0, indicative)")
print(f"candidates: {len(cands)} in {nframes} frames "
      f"({per_frame.mean():.1f}/frame, min {per_frame.min()}, max {per_frame.max()})")
print(f"ROIs      : {len(rois)} kept, {len(cands) - len(rois)} dropped at the border")
if len(rois):
    print(f"ROI photons: median sum={np.median(rois.images.sum(axis=(1, 2))):.0f}")
    # the emitter should sit at the ROI centre: check the centre-of-mass
    r = rois.images - np.median(photons)
    yy, xx = np.mgrid[0:roisize, 0:roisize]
    w = np.clip(r, 0, None).sum(axis=(1, 2))
    cx = (np.clip(r, 0, None) * xx).sum(axis=(1, 2)) / w
    cy = (np.clip(r, 0, None) * yy).sum(axis=(1, 2)) / w
    print(f"ROI centre of mass: x={np.median(cx):.2f} y={np.median(cy):.2f} "
          f"(expected {(roisize - 1) / 2:.1f})")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("\n(matplotlib not available: no figure written)")
    raise SystemExit(0)

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
f0 = cands.frame == start
ax[0].imshow(photons[0], cmap="gray"); ax[0].set_title("photons, frame 0")
ax[1].imshow(filtered[0], cmap="gray")
ax[1].plot(cands.x[f0], cands.y[f0], "o", mfc="none", mec="r", ms=10)
ax[1].set_title(f"filtered + {f0.sum()} candidates")
n = min(64, len(rois))
if n:
    side = int(np.ceil(np.sqrt(n)))
    montage = np.zeros((side * roisize, side * roisize), np.float32)
    for i in range(n):
        r, c = divmod(i, side)
        montage[r*roisize:(r+1)*roisize, c*roisize:(c+1)*roisize] = rois.images[i]
    ax[2].imshow(montage, cmap="gray"); ax[2].set_title(f"first {n} ROIs")
for a in ax: a.set_xticks([]); a.set_yticks([])
out = Path("detection_check.png"); fig.tight_layout(); fig.savefig(out, dpi=110)
print(f"\nwrote {out.resolve()}")
