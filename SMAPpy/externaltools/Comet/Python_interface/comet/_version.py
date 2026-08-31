"""Single source of truth for the package version.

Kept in its own dependency-free module so setuptools can read it statically at
build time (see ``[tool.setuptools.dynamic]`` in pyproject.toml) without
importing the rest of the package.
"""

__version__ = "1.1.0"
