"""Lazy wrappers around the Tk file dialogs.

tkinter is part of the standard library but is not present in every Python
installation: many Linux distributions ship it as a separate ``python3-tk``
package and most slim container images leave it out entirely.  Importing it at
module level would make ``import comet`` fail on those systems even though the
dialogs are only ever needed when a path is omitted in an interactive session,
so it is imported on demand instead.
"""

_MISSING_TK_MSG = (
    "No file name was given, so COMET tried to open an interactive file "
    "dialog, but tkinter is not available in this Python installation. "
    "Pass an explicit path instead, or install tkinter "
    "(e.g. 'sudo apt install python3-tk' on Debian/Ubuntu)."
)


def _filedialog():
    """Return the tkinter.filedialog module with a hidden root window."""
    try:
        import tkinter
        from tkinter import filedialog
    except ImportError as exc:  # pragma: no cover - depends on the interpreter
        raise RuntimeError(_MISSING_TK_MSG) from exc
    tkinter.Tk().withdraw()
    return filedialog


def ask_open_filename(**kwargs):
    """Prompt for an existing file. Raises RuntimeError if tkinter is missing."""
    return _filedialog().askopenfilename(**kwargs)


def ask_save_filename(**kwargs):
    """Prompt for a save location. Raises RuntimeError if tkinter is missing."""
    return _filedialog().asksaveasfilename(**kwargs)


def ask_directory(**kwargs):
    """Prompt for a directory. Raises RuntimeError if tkinter is missing."""
    return _filedialog().askdirectory(**kwargs)
