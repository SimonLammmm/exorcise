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

__version__ = "3.0.0"

__all__ = ["__version__"]
