"""Small shared helpers.

These were previously scattered between crispr_tools.tools (a 722-line module
that also pulled in scipy, h5py, statsmodels and xlsxwriter at import time) and
crispr_screen_viewer.functions_etc (which imported Dash). Collecting the dozen
functions that are actually used means neither of those heavyweight imports
happens any more.
"""

from __future__ import annotations

import os
import sys
import unicodedata
from pathlib import Path
from typing import Dict, Iterator, Union

import pandas as pd
from loguru import logger

# Joins the two halves of a comparison in output filenames and column names.
ARROW = "→"

# Human-readable names for the timepoint conventions used in screen metadata.
TIMEPOINT_LABELS = {
    "fromstart": "From experiment start",
    "otherprior": "From midpoint",
    "endpoints": "Matched time points",
}

# Which column of an analysis result table carries the effect size.
ANALYSIS_SCORE_NAMES = {"drugz": "normZ", "mageck": "lfc", "chronos": "chronos_score"}

STATISTIC_LABELS = {
    "lfc": "Log2(FC)",
    "normZ": "NormZ",
    "chronos_score": "Chronos score",
    "neg_p": "Sensitising p-value",
    "pos_p": "Suppressing p-value",
    "neg_fdr": "Sensitising FDR",
    "pos_fdr": "Suppressing FDR",
    "fdr": "FDR",
    "p": "p-value",
    "fdr_log10": "-log10(FDR)",
}


#### Paths and files ####


def maybe_its_gz(filename) -> str:
    """Return "<filename>.gz" when filename is absent but the gzip is present."""
    filename = str(filename)
    if os.path.isfile(filename):
        return filename
    gzipped = filename + ".gz"
    if os.path.isfile(gzipped):
        return gzipped
    return filename


def file_exists_and_not_empty(fn) -> bool:
    """Checkpoint test: True when a previous run already produced this file."""
    return os.path.exists(fn) and os.path.getsize(fn) != 0


def is_temp_file(fn: Union[str, Path]) -> bool:
    """True for editor backups and the dotfiles that macOS scatters over shares."""
    fn = os.path.basename(str(fn))
    return fn.startswith(".") or fn.startswith("~")


def iter_files_by_prefix(prefix: Union[str, Path], req_suffix=None,
                         allow_gz=True) -> Iterator[str]:
    """Yield the names of files in prefix's directory whose names contain the
    prefix's final component, optionally filtered by suffix."""
    prefix = Path(prefix)

    if req_suffix is None:
        def matches_suffix(name: str) -> bool:
            return True
    elif allow_gz:
        def matches_suffix(name: str) -> bool:
            return name.endswith(req_suffix) or name.endswith(req_suffix + ".gz")
    else:
        def matches_suffix(name: str) -> bool:
            return name.endswith(req_suffix)

    for name in os.listdir(prefix.parent):
        if matches_suffix(name) and prefix.parts[-1] in name and not is_temp_file(name):
            yield name


def ensure_parent_dir(path: Union[str, Path]) -> None:
    """Create the directory a file is about to be written into.

    os.path.dirname("out") is "", which makes a bare os.makedirs() raise. That
    caught out anyone who passed a --prefix without a directory component.
    """
    parent = os.path.dirname(os.path.expanduser(str(path)))
    if parent:
        os.makedirs(parent, exist_ok=True)


#### Values ####


def list_not_str(thing):
    """Wrap a bare string in a list, leave other iterables alone."""
    if isinstance(thing, str):
        return [thing]
    return thing


def is_nt(s) -> bool:
    """True when a treatment field means "no treatment"."""
    return pd.isna(s) or (s in {"None", "DMSO", "", "WT", "NT"})


def normalise_text(s: str) -> str:
    """Fold accents and formatting characters out of free text."""
    return unicodedata.normalize("NFKD", s)


def clargs_to_dict(argstring: str) -> Dict[str, str]:
    """Parse "--arg1 val --flag --arg2 val2" into {'arg1': 'val', 'flag': '',
    'arg2': 'val2'}.

    Used for the Arguments column of an analysis workbook. It lives here rather
    than in pipeline.py because the workbook parser needs it too, and importing
    it from there created a circular import between the two modules.
    """
    kwargs: Dict[str, str] = {}
    parts = argstring.split()
    for i, part in enumerate(parts):
        if not part.startswith("-"):
            continue
        key = part.lstrip("-")
        if i + 1 == len(parts) or parts[i + 1].startswith("-"):
            kwargs[key] = ""
        else:
            kwargs[key] = parts[i + 1]
    return kwargs


#### Tables ####


def df_rename_columns(df: pd.DataFrame, newcols: dict, inplace=False,
                      axis="columns") -> pd.Index:
    """Rename only the labels present in newcols, leaving the rest as they are."""
    if axis in ("columns", 0):
        axis = "columns"
    elif axis in ("index", 1):
        axis = "index"
    else:
        raise ValueError('axis needs to be "columns"|0 or "index"|1')

    mapper = {k: k for k in getattr(df, axis)}
    mapper.update(newcols)
    renamed = getattr(df, axis).map(mapper)

    if not inplace:
        return renamed
    setattr(df, axis, renamed)
    return renamed


def load_stats_csv(
    fn,
    drop_controls_with=("control", "Control", "CONTROL", "Non-target", "Cutting"),
) -> pd.DataFrame:
    """Load a double-headered statistics CSV, dropping control rows."""
    df = pd.read_csv(fn, index_col=0, header=[0, 1])

    if df.index.isna().any():
        df = df.loc[~df.index.isna()]
        logger.warning(f"NaN genes found in (and dropped from) {fn}")

    if drop_controls_with:
        drop = df.index != df.index  # all False, right length
        for indicator in drop_controls_with:
            drop = drop | df.index.str.contains(indicator)
        df = df.loc[~drop]

    return df


#### Logging ####


def set_loguru_level(target_logger, level="INFO") -> None:
    """Replace loguru's default sink with one at the requested level."""
    try:
        target_logger.remove(0)
    except ValueError:
        pass
    target_logger.add(sys.stderr, level=level)
