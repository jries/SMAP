import subprocess
import sys
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

extra_compile_args = ["-O3"]
extra_link_args = []


def _macos_libcxx_workaround():
    """Work around a Command Line Tools install with incomplete libc++ headers.

    Some macOS setups have an almost empty ``CommandLineTools/usr/include/c++/v1``
    while the SDK's copy is complete, so any ``#include <cstddef>`` fails.  When
    that is the case, point the compiler at the SDK's headers instead.

    The real fix is to repair the toolchain, e.g. by pointing xcode-select at a
    full Xcode (``sudo xcode-select -s /Applications/Xcode.app/Contents/Developer``)
    or reinstalling the Command Line Tools.
    """
    if sys.platform != "darwin":
        return []
    if _compiles([]):
        return []  # the toolchain is fine, which is the normal case

    try:
        developer = subprocess.check_output(["xcode-select", "-p"], text=True).strip()
    except Exception:
        developer = ""

    candidates = [Path("/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")]
    if developer:
        candidates.insert(0, Path(developer) / "SDKs/MacOSX.sdk")
    for sdk in candidates:
        headers = sdk / "usr/include/c++/v1"
        flags = ["-nostdinc++", "-isystem", str(headers)]
        if (headers / "cstddef").exists() and _compiles(flags):
            print(f"note: incomplete libc++ in the active toolchain; "
                  f"using the headers from {headers}")
            return flags

    print("warning: no working C++ standard library found; the build will "
          "probably fail. Try 'sudo xcode-select --install', or point "
          "xcode-select at a full Xcode and run 'sudo xcodebuild -license accept'.")
    return []


def _compiles(flags):
    """Whether a trivial C++ program compiles with these extra flags."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        source = Path(tmp) / "probe.cpp"
        source.write_text("#include <cstddef>\n#include <thread>\nint main(){}\n")
        try:
            return subprocess.run(
                ["c++", *flags, str(source), "-o", str(Path(tmp) / "probe")],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            ).returncode == 0
        except Exception:
            return False


extra_compile_args += _macos_libcxx_workaround()

ext_modules = [
    Pybind11Extension(
        "smapfit._fit3d",
        ["csrc/fit.cpp"],
        include_dirs=["csrc"],
        cxx_std=17,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
