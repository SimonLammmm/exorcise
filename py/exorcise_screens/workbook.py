"""Parse the Excel workbook that describes how a screen should be analysed.

Was the AnalysisWorkbook class in crispr_tools/data_classes.py. Three of its six
methods were unreachable (verify_wb_integrity, iter_comps, check_analyses_ran)
and have been dropped along with the rest of that module.

The workbook has five sheets:

  Experiment details   two columns, field name and value
  Sample details       one row per replicate
  Control groups       which sample is the control for which
  Analyses             one row per analysis to run
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from loguru import logger

from exorcise_screens.util import clargs_to_dict

# Mapping from the labels used in the workbook to the internal keys the pipeline
# expects. Several spellings are accepted because workbooks in the wild use the
# internal style, and older ones use different labels again.
EXPERIMENT_DETAIL_KEYS = {
    "Experiment name": "experiment_id",
    "Analysis version": "analysis_version",
    "Notes": "notes",
    "File prefix": "file_prefix",
    # Some workbooks are written in the internal style.
    "analysis_version": "analysis_version",
    "file_prefix": "file_prefix",
    # Legacy labels.
    "Analysis name": "experiment_id",
    "Library": "library",
}

REQUIRED_SHEETS = (
    "Experiment details",
    "Sample details",
    "Control groups",
    "Analyses",
)


def _safesplit(value: str) -> List[str]:
    """Split a comma-separated cell, tolerating spaces around the commas."""
    return str(value).replace(" ", "").split(",")


def _fill_blanks(frame):
    """Replace NaN with the empty string, keeping every column as object.

    Sheets are read with dtype=object on purpose: nothing in a workbook needs to
    be numeric, and letting pandas guess turns sample names like "01" into 1.
    Plain fillna("") would silently downcast the result and, since pandas 2.2,
    warn about doing so on every read.
    """
    try:
        with pd.option_context("future.no_silent_downcasting", True):
            return frame.fillna("")
    except pd.errors.OptionError:
        # pandas older than 2.2 has no such option and does not downcast here.
        return frame.fillna("")


def read_experiment_id(fn) -> str:
    """The experiment's name, read from the Experiment details sheet alone.

    Deliberately lighter than constructing an AnalysisWorkbook. Removing an
    experiment from a database needs only its ID, and should not fail because
    some sheet it will never look at is malformed or absent.
    """
    details = _fill_blanks(
        pd.read_excel(fn, sheet_name="Experiment details", header=None, index_col=0)[1]
    )

    # "Analysis name" is the legacy label for the same field.
    for label in ("Experiment name", "Analysis name"):
        if label in details.index:
            value = str(details[label]).strip()
            if value:
                return value

    raise KeyError(
        f'{fn} has no non-empty "Experiment name" row in its "Experiment '
        f'details" sheet, so the experiment cannot be identified.'
    )


class AnalysisWorkbook:
    """An analysis workbook, and the configuration dictionary it describes."""

    def __init__(self, xlfn, parse_workbook=True, counts_dir="."):
        self.fn = xlfn
        self.wb = self.load_analysis_workbook(xlfn)

        self.experiment_details: pd.Series = self.wb["Experiment details"]
        self.samples: pd.DataFrame = self.wb["Sample details"]
        self.replicates: pd.DataFrame = self.wb["Replicate details"]
        self.control_groups: pd.DataFrame = self.wb["Control groups"]
        self.analyses: pd.DataFrame = self.wb["Analyses"]

        self.counts_dir = counts_dir
        self.expd: Dict[str, Any] = {}
        if parse_workbook:
            self.expd = self.workbook_to_dict(counts_dir=counts_dir)

    # Retained because update_database and older callers spell it this way.
    @property
    def analylses(self) -> pd.DataFrame:  # noqa: N802 - historical misspelling
        """Deprecated misspelling of `analyses`."""
        return self.analyses

    @staticmethod
    def load_analysis_workbook(fn) -> Dict[str, Any]:
        """Read every sheet. Nothing is coerced to a number at this stage."""
        wb = pd.read_excel(fn, sheet_name=None, dtype=object)
        wb = {name: _fill_blanks(sheet) for name, sheet in wb.items()}

        missing = [s for s in REQUIRED_SHEETS if s not in wb]
        if missing:
            raise KeyError(
                f"{fn} is missing the sheet(s): {', '.join(missing)}. "
                f"An analysis workbook needs: {', '.join(REQUIRED_SHEETS)}."
            )

        # "Sample details" is per replicate on disk. Keep both views: one row per
        # replicate, and one row per sample.
        wb["Replicate details"] = wb["Sample details"].set_index("Replicate", drop=False)
        wb["Sample details"] = (
            wb["Replicate details"].drop_duplicates("Sample").set_index("Sample", drop=False)
        )

        # Experiment details is a two-column key/value sheet, not a table.
        wb["Experiment details"] = _fill_blanks(
            pd.read_excel(fn, sheet_name="Experiment details", header=None,
                          index_col=0)[1]
        )

        return wb

    def workbook_to_dict(self, counts_dir="") -> Dict[str, Any]:
        """Build the configuration dictionary the pipeline runs from."""
        counts_dir = Path(counts_dir)
        logger.info(f"Parsing workbook {self.fn}")

        config: Dict[str, Any] = {}

        details = self.wb["Experiment details"].to_dict()
        for label, key in EXPERIMENT_DETAIL_KEYS.items():
            if label in details:
                config[key] = details[label]

        config["control_groups"] = self._control_groups()
        config.update(self._sample_groupings())
        config["analyses"] = self._analyses(counts_dir)

        if "experiment_id" not in config:
            raise KeyError(
                f'{self.fn} "Experiment details" has no "Experiment name" row, '
                "which is what identifies the experiment."
            )
        config.setdefault("file_prefix", config["experiment_id"])

        return config

    def _control_groups(self) -> Dict[str, Dict[str, list]]:
        """{group: {control_sample: [test_sample, ...]}}"""
        groups = self.wb["Control groups"]
        result = {}
        for name in groups.Group.unique():
            in_group = groups[groups.Group == name]
            result[name] = {
                control: list(in_group.loc[rows, "Test sample"].values)
                for control, rows in in_group.groupby("Control sample").groups.items()
            }
        return result

    def _sample_groupings(self) -> Dict[str, Dict[str, list]]:
        """Replicate, timepoint and trajectory groupings, keyed by sample."""
        replicates = self.wb["Replicate details"]

        groupings = {
            "sample_reps": {
                k: list(v.values)
                for k, v in replicates.groupby("Sample").groups.items()
            }
        }

        # Chronos needs to know how long each replicate grew, and which
        # replicates are the same sample at different timepoints.
        for column, key in (("Days grown", "days_grown"), ("Trajectory", "cell_line_hash")):
            if column not in replicates.columns:
                logger.info(
                    f'"Sample details" has no "{column}" column. '
                    "Chronos analyses will not be possible."
                )
                groupings[key] = {}
                continue
            indexed = replicates.set_index(replicates[column].astype(str))
            groupings[key] = {
                k: list(v.values) for k, v in indexed.groupby("Sample").groups.items()
            }

        return groupings

    def _analyses(self, counts_dir: Path) -> List[Dict[str, Any]]:
        """One dictionary per analysis to run.

        A single row can name several methods, in which case it becomes one
        analysis per method.
        """
        analyses: List[Dict[str, Any]] = []

        for _, row in self.wb["Analyses"].iterrows():
            for method in _safesplit(row["Method"]):
                pseudocount = int(row["Add pseudocount"]) if row["Add pseudocount"] else 1

                analysis: Dict[str, Any] = {
                    "method": method,
                    "groups": _safesplit(row["Control group"]),
                    "pseudocount": pseudocount,
                    "kwargs": {},
                    "name": row["Name"] if row.get("Name") else "-",
                }

                if row["Counts file"]:
                    analysis["counts_file"] = os.path.join(counts_dir, row["Counts file"])

                analysis["kwargs"] = self._analysis_kwargs(row, method, pseudocount)
                self._apply_paired(analysis, row, method)

                analyses.append(analysis)

        return analyses

    @staticmethod
    def _analysis_kwargs(row, method: str, pseudocount: int) -> Dict[str, Any]:
        """The Arguments cell, as keyword arguments for the analysis program."""
        arguments = row["Arguments"]

        if not arguments:
            # MAGeCK's default: discard guides that are absent from both arms.
            if method == "mageck":
                cutoff = 2 * pseudocount if pseudocount > 1 else 0
                return {"remove-zero": "both", "remove-zero-threshold": cutoff}
            return {}

        # The cell may hold a Python dict literal or a command line fragment.
        try:
            kwargs = eval(arguments)  # noqa: S307 - workbooks are trusted input
            if not isinstance(kwargs, dict):
                raise ValueError("not a dict")
        except Exception:
            logger.info(f"Parsing Arguments {arguments!r} as command line options")
            kwargs = clargs_to_dict(str(arguments))

        logger.info(f"Arguments {arguments!r} parsed to {kwargs}")
        return kwargs

    @staticmethod
    def _apply_paired(analysis: Dict[str, Any], row, method: str) -> None:
        """MAGeCK and DrugZ spell paired analysis in opposite directions."""
        paired = row["Paired"]
        if pd.isna(paired):
            return
        if paired and method == "mageck":
            analysis["kwargs"]["paired"] = ""
        elif not paired and method == "drugz":
            analysis["kwargs"]["unpaired"] = True
