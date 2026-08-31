"""Runtime detection of the available COMET compute backends.

COMET can run on three backends:

``cuda``
    numba-cuda kernels; needs an NVIDIA GPU and a working CUDA driver.
``torch``
    PyTorch kernels; uses CUDA or Apple MPS when present, otherwise runs on
    the CPU. Requires the optional ``torch`` dependency.
``cpu``
    Pure NumPy fallback; always available, but slow on large datasets.

Everything here is defensive: probing a GPU driver can fail in a lot of
creative ways, and a detection helper must never be the thing that raises.
"""

BACKENDS = ("cuda", "torch", "cpu")


def cuda_available():
    """True if numba can talk to a CUDA device."""
    try:
        from numba import cuda
        return bool(cuda.is_available())
    except Exception:
        # Missing driver, missing toolkit, unsupported arch, ...
        return False


def torch_available():
    """True if PyTorch is importable."""
    try:
        import torch  # noqa: F401
        return True
    except Exception:
        return False


def torch_device():
    """Return the best torch device name, or None if torch is unusable."""
    try:
        import torch
    except Exception:
        return None
    try:
        if torch.cuda.is_available():
            return "cuda"
    except Exception:
        pass
    try:
        # torch.backends.mps only exists on torch >= 1.12
        if torch.backends.mps.is_available():
            return "mps"
    except Exception:
        pass
    return "cpu"


def torch_accelerated():
    """True if torch is installed and backed by a GPU worth preferring.

    Deliberately CUDA-only. Apple MPS is a real GPU -- it sustains ~2.8 TFLOP/s
    on matmul here, comfortably beating the CPU -- but it loses badly on *this*
    workload. One cost+gradient evaluation, M-series machine, across the range
    of segment counts real acquisitions produce:

        pairs        segments    numba cpu    torch mps
          496 k           150      3.4 ms       23.8 ms
        4 048 k         1 000     19.7 ms      108.6 ms
       12 545 k         2 500     57.0 ms      349.5 ms

    The reason is that eager-mode torch materialises every intermediate across
    all pairs. At 4 M pairs each (P, 3) temporary is ~48 MB, and there are a
    dozen of them, so the evaluation is bound by memory traffic rather than
    arithmetic: ~55% of the time is gathers and elementwise maths, ~45% the
    index_add_ that accumulates the gradient. The compiled CPU kernel keeps each
    pair in registers and never materialises any of it.

    Competing would need a *fused* GPU kernel, which is exactly what the
    numba-cuda backend already is for NVIDIA hardware.

    MPS is also float32-only, but that is not the reason: per-evaluation
    gradients differ by up to ~1e-5 at 3 M pairs while the converged drift
    agrees with the CPU backend to 0.0013 nm. This is a speed decision.

    MPS stays selectable with mode="torch"; it is just not the automatic choice.
    """
    return torch_device() == "cuda"


def available_backends():
    """List the backends usable on this machine, fastest first."""
    backends = []
    if cuda_available():
        backends.append("cuda")
    if torch_available():
        backends.append("torch")
    backends.append("cpu")
    return backends


def best_backend():
    """Pick the fastest usable backend.

    Order: numba-cuda, then torch if it has CUDA, otherwise the compiled CPU
    kernel. Neither a CPU-only nor an MPS torch install is preferred over the
    CPU backend -- both are slower than it (see :func:`torch_accelerated`).
    """
    if cuda_available():
        return "cuda"
    if torch_accelerated():
        return "torch"
    return "cpu"


def describe_backends():
    """Human-readable one-liner per backend, for diagnostics."""
    lines = []
    lines.append("cuda  : {}".format(
        "available" if cuda_available() else "unavailable (no NVIDIA GPU / CUDA driver)"))
    if torch_available():
        device = torch_device()
        note = ""
        if device == "mps":
            note = "; slower than cpu for this workload, not auto-selected"
        elif device == "cpu":
            note = "; no GPU, not auto-selected"
        lines.append("torch : available (device: {}{})".format(device, note))
    else:
        lines.append("torch : unavailable (pip install py-comet[torch])")
    lines.append("cpu   : available (numba-compiled)")
    return lines
