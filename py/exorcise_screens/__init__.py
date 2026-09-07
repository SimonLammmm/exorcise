"""Screen analysis tools that ship alongside Exorcise.

This package is the merger of two projects:

* crispr_tools (John C. Thomas), which provided read counting and the
  MAGeCK/DrugZ/Chronos analysis pipeline;
* crispr_screen_viewer (John C. Thomas), whose database builder is kept here.
  Its Dash web application is not part of this package.

Only the functionality reachable from the three command line entry points is
carried over. See CHANGELOG.md for what was dropped and why.

Nothing is imported eagerly here: the analysis pipeline pulls in TensorFlow by
way of Chronos, and read counting should not have to wait for that.
"""

def _read_version() -> str:
    """The project version.

    The single source of truth is the VERSION file at the root of the
    installation, which is also what the R entry points and the image tag read.

    That file is consulted first, because pyproject.toml derives the package
    version from this attribute: at build time the package metadata does not
    exist yet, so the file is the only thing available. Once installed the file
    is usually absent, and the metadata written from it at build time is used
    instead.
    """
    from pathlib import Path

    # py/exorcise_screens/__init__.py -> py/exorcise_screens -> py -> root
    version_file = Path(__file__).resolve().parents[2] / "VERSION"
    try:
        version = version_file.read_text(encoding="utf-8").strip()
        if version:
            return version
    except OSError:
        pass

    try:
        from importlib.metadata import PackageNotFoundError, version as _metadata

        return _metadata("exorcise-screens")
    except PackageNotFoundError:
        return "unknown"


__version__ = _read_version()

__all__ = ["__version__"]
