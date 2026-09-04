"""Run batches of CRISPR screen analyses described by a configuration file.

Three methods are available: MAGeCK (an external command), DrugZ (called in
process, vendored in this package) and Chronos (called in process, needs
TensorFlow). Which of them run is decided entirely by the Analyses sheet of the
workbook, or the `analyses` key of a JSON configuration.

Every analysis is checkpointed on its output file, so a rerun only does the work
that is missing unless --overwrite is given.

Was crispr_tools/crispr_pipeline.py.
"""

from __future__ import annotations

import datetime
import gzip
import json
import math
import os
import shlex
import subprocess
from copy import copy
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from attrdictionary import AttrDict
from loguru import logger

from exorcise_screens.drugz import drugZ_analysis
from exorcise_screens.tables import TABULATORS
from exorcise_screens.util import (
    file_exists_and_not_empty,
    list_not_str,
    maybe_its_gz,
)
from exorcise_screens.workbook import AnalysisWorkbook

#: Joins control and treatment sample names in output filenames.
FILENAME_JOINER = "-"

#: Guides matching these are treated as negative controls by Chronos. The first
#: pattern excludes Exorcise's own exo_Non-targeting labels so that the library
#: designer's declared controls take precedence.
CHRONOS_DESIGN_CONTROLS = r".*(?<!exo_)Non-targeting.*|.*_OR\d.*|.*_Cutting\d.*"
CHRONOS_FALLBACK_CONTROLS = r".*Non-targeting.*|.*_OR\d.*|.*_Cutting\d.*"


class PipelineOptionsError(Exception):
    """A configuration file problem that stops the pipeline before it starts."""


#### DrugZ ####


def call_drugz_batch(sample_reps: Dict[str, list], days_grown, cell_line_hash,
                     control_map: Dict[str, list], counts_file: str, prefix: str,
                     kwargs: Optional[dict] = None, pseudocount=1,
                     drop_guide_less_than=0, drop_guide_type=None,
                     overwrite=False) -> None:
    """One DrugZ result file per comparison, at {prefix}.{ctrl}-{treat}.tsv."""
    kwargs = dict(kwargs or {})

    # These two can arrive either as arguments or through the workbook's
    # Arguments column.
    if "drop_guide_less_than" in kwargs:
        drop_guide_less_than = int(kwargs.pop("drop_guide_less_than"))
    if "drop_guide_type" in kwargs:
        drop_guide_type = kwargs.pop("drop_guide_type")

    if drop_guide_less_than and drop_guide_type not in ("any", "all", "both"):
        raise PipelineOptionsError(
            f"drop_guide_type must be 'any', 'all' or 'both', got {drop_guide_type!r}."
        )

    args = AttrDict()
    args.infile = counts_file
    args.index_column = 0
    args.minobs = 1
    args.half_window_size = 500
    args.quiet = False
    args.pseudocount = pseudocount
    args.fc_outfile = ""
    args.unpaired = False
    args.remove_genes = None
    args.update(kwargs)

    # Dropping low-abundance guides means writing a reduced counts file per
    # comparison, since DrugZ reads from disk.
    temp_counts = None
    counts = None
    if drop_guide_less_than:
        temp_counts = f"{prefix}.tmp_drugz_counts.tsv"
        args.infile = temp_counts
        counts = pd.read_csv(counts_file, index_col=0, sep="\t")

    try:
        for control_sample, treat_samples in control_map.items():
            args.control_samples = ",".join(list_not_str(sample_reps[control_sample]))

            for treat_sample in list_not_str(treat_samples):
                args.drug_samples = ",".join(list_not_str(sample_reps[treat_sample]))
                args.drugz_output_file = os.path.normpath(
                    f"{prefix}.{control_sample}{FILENAME_JOINER}{treat_sample}.tsv"
                )

                if not overwrite and file_exists_and_not_empty(
                    maybe_its_gz(args.drugz_output_file)
                ):
                    logger.info(
                        f"DrugZ result exists, skipping: {args.drugz_output_file}"
                    )
                    continue

                if drop_guide_less_than:
                    _write_filtered_counts(
                        counts, temp_counts, sample_reps, control_sample,
                        treat_sample, drop_guide_less_than, drop_guide_type,
                    )

                drugZ_analysis(args)
    finally:
        if temp_counts and os.path.exists(temp_counts):
            os.remove(temp_counts)

    logger.info("Finished DrugZ")


def _write_filtered_counts(counts: pd.DataFrame, temp_fn: str, sample_reps, control,
                           treat, threshold: int, mode: str) -> None:
    """Write a counts file with sparsely-observed guides removed."""
    replicates = list_not_str(sample_reps[control]) + list_not_str(sample_reps[treat])
    below = counts.loc[:, replicates] < threshold
    drop = below.any(axis=1) if mode == "any" else below.all(axis=1)

    counts.loc[~drop, ["gene"] + replicates].to_csv(temp_fn, sep="\t")
    logger.info(
        f"Dropped {int(drop.sum())} guides from {control}{FILENAME_JOINER}{treat} "
        f"scoring under {threshold} in '{mode}' replicate(s)."
    )


#### MAGeCK ####


def mageck_version() -> str:
    """Confirm MAGeCK is installed and report its version."""
    try:
        return subprocess.check_output(["mageck", "-v"], text=True).strip()
    except FileNotFoundError:
        raise PipelineOptionsError(
            "MAGeCK is not on the PATH. Install it, or drop 'mageck' from the "
            "Method column of your workbook."
        ) from None


def call_mageck(control_sample: str, treat_sample: str,
                sample_reps: Dict[str, List[str]], counts_file: str, prefix: str,
                kwargs: Optional[Dict] = None, dryrun=False) -> None:
    """Run one MAGeCK comparison.

    The command is assembled as an argument list rather than a shell string.
    Previously it was joined into one string and run with shell=True, which meant
    a sample name or path containing a space or a shell metacharacter would break
    the command or inject into it.
    """
    command = [
        "mageck", "test",
        "-k", str(counts_file),
        "-t", ",".join(sample_reps[treat_sample]),
        "-c", ",".join(sample_reps[control_sample]),
        "-n", f"{prefix}.{control_sample}{FILENAME_JOINER}{treat_sample}",
    ]

    for key, value in (kwargs or {}).items():
        command.append(f"--{key}")
        # An empty value means a bare flag, e.g. {'paired': ''}.
        if value != "" and value is not None:
            command.append(str(value))

    logger.info("Running: " + " ".join(shlex.quote(part) for part in command))
    if dryrun:
        return

    result = subprocess.run(command)
    if result.returncode != 0:
        logger.error(
            f"MAGeCK exited with status {result.returncode} for "
            f"{control_sample}{FILENAME_JOINER}{treat_sample}."
        )


def call_mageck_batch(sample_reps: Dict[str, list], days_grown, cell_line_hash,
                      control_map: Dict[str, list], counts_file: str, prefix: str,
                      kwargs: Optional[dict] = None, pseudocount=1,
                      overwrite=False) -> None:
    """Run every comparison in control_map through MAGeCK."""
    logger.info(f"Running MAGeCK version {mageck_version()}")

    # MAGeCK cannot read gzip and has no pseudocount option, so both are handled
    # by writing a temporary counts file.
    temp_counts = prefix + ".tmp_mageck_counts.tsv"
    counts = pd.read_csv(maybe_its_gz(counts_file), sep="\t", index_col=0)
    if pseudocount > 1:
        numeric = counts.dtypes != object
        counts.loc[:, numeric] += pseudocount
    counts.to_csv(temp_counts, sep="\t")

    try:
        for control_sample, treat_samples in control_map.items():
            for treat_sample in list_not_str(treat_samples):
                if treat_sample == control_sample:
                    continue

                outfn = (
                    f"{prefix}.{control_sample}{FILENAME_JOINER}{treat_sample}"
                    ".gene_summary.txt"
                )
                if not overwrite and file_exists_and_not_empty(maybe_its_gz(outfn)):
                    logger.info(f"MAGeCK result exists, skipping: {outfn}")
                    continue

                call_mageck(
                    control_sample, treat_sample, sample_reps, temp_counts,
                    prefix, kwargs,
                )
    finally:
        if os.path.exists(temp_counts):
            os.remove(temp_counts)


#### Chronos ####


def call_chronos_batch(sample_reps: Dict[str, list], days_grown: Dict[str, list],
                       cell_line_hash: Dict[str, list], control_map: Dict[str, list],
                       counts_file: str, prefix: str, kwargs: Optional[dict] = None,
                       pseudocount=1, overwrite=False) -> None:
    """Fit a Chronos model per control group.

    Chronos models growth over time, so it only applies to comparisons spanning
    more than one timepoint where the control is the earliest one.
    """
    counts = pd.read_csv(maybe_its_gz(counts_file), sep="\t")
    if pseudocount > 1:
        numeric = counts.dtypes != object
        counts.loc[:, numeric] += pseudocount

    guide_column, gene_column = counts.columns[0], counts.columns[1]

    readcounts = (
        counts.drop(gene_column, axis=1)
        .set_index(guide_column)
        .transpose()
        .reset_index()
        .rename(columns={"index": "sequence_ID"})
    )

    guidemap = pd.DataFrame({"sgrna": counts[guide_column], "gene": counts[gene_column]})
    negative_controls = _chronos_negative_controls(counts[guide_column])
    sequence_map = _chronos_sequence_map(sample_reps, days_grown, cell_line_hash)

    for control_sample, treat_samples in control_map.items():
        group = [control_sample] + list(list_not_str(treat_samples))
        group_map = sequence_map[sequence_map["sample"].isin(group)].copy()
        group_map.loc[group_map["sample"] == control_sample, "cell_line_name"] = "pDNA"
        group_map = group_map.drop_duplicates(subset=["sequence_ID"])

        timepoints = set(group_map["days"])
        if len(timepoints) == 1:
            logger.info(
                f"Skipping Chronos for {group}: every sample is from the same "
                "timepoint."
            )
            continue

        control_timepoints = set(
            group_map.loc[group_map["cell_line_name"] == "pDNA", "days"]
        )
        if not control_timepoints or min(timepoints) != min(control_timepoints):
            logger.info(
                f"Skipping Chronos for {group}: the control sample is not from "
                "the earliest timepoint."
            )
            continue

        group_prefix = f"{prefix}.{control_sample}"
        if not overwrite and file_exists_and_not_empty(
            maybe_its_gz(os.path.join(group_prefix, "gene_effect.hdf5"))
        ):
            logger.info(f"Chronos result exists, skipping: {group_prefix}")
            continue

        _fit_chronos(
            group_map, readcounts, guidemap, negative_controls, group_prefix,
            min(timepoints), control_sample,
        )


def _chronos_negative_controls(guides: pd.Series) -> pd.Series:
    """Prefer the library designer's controls; fall back to Exorcise's labels."""
    designed = guides[guides.str.match(CHRONOS_DESIGN_CONTROLS)]
    if len(designed):
        return designed
    return guides[guides.str.match(CHRONOS_FALLBACK_CONTROLS)]


def _unstack_to_frame(mapping: Dict[str, list], value_name: str) -> pd.DataFrame:
    """Turn {sample: [values]} into a two-column sample/value frame."""
    frame = (
        pd.DataFrame(pd.DataFrame.from_dict(mapping, orient="index").unstack())
        .reset_index()
        .drop("level_0", axis=1)
        .rename(columns={"level_1": "sample", 0: value_name})
    )
    return frame.loc[[v is not None for v in frame[value_name]]]


def _chronos_sequence_map(sample_reps, days_grown, cell_line_hash) -> pd.DataFrame:
    """The replicate/timepoint/cell-line table Chronos wants."""
    sequences = _unstack_to_frame(sample_reps, "sequence_ID")

    days = _unstack_to_frame(days_grown, "days")
    days["days"] = days["days"].astype(float)
    days = days.loc[[not math.isnan(d) for d in days["days"]]]

    lines = _unstack_to_frame(cell_line_hash, "cell_line_name")
    lines["cell_line_name"] = lines["cell_line_name"].astype(str)

    sequence_map = (
        sequences.merge(days, on="sample", how="left")
        .drop_duplicates(subset=["sample"])
        .merge(lines, on="sample", how="left")
        .drop_duplicates(subset=["sample"])
    )
    sequence_map["pDNA_batch"] = 0
    return sequence_map


def _fit_chronos(group_map, readcounts, guidemap, negative_controls, group_prefix,
                 initial_screen_delay, control_sample) -> None:
    """Train and save one Chronos model.

    Imported lazily: Chronos pulls in TensorFlow, and a MAGeCK-only run should
    not pay for that.
    """
    import chronos

    # A cell line with no replicate after the control timepoint contributes
    # nothing and upsets the fit.
    early = group_map.loc[group_map["days"] <= initial_screen_delay, "cell_line_name"]
    for line in early[early != "pDNA"]:
        line_days = group_map.loc[group_map["cell_line_name"] == line, "days"]
        if not any(day > initial_screen_delay for day in line_days):
            group_map = group_map.loc[group_map["cell_line_name"] != line]

    group_counts = readcounts[
        readcounts["sequence_ID"].isin(group_map["sequence_ID"])
    ].set_index("sequence_ID")
    chronos.nan_outgrowths(group_counts, group_map, guidemap)

    common = dict(
        readcounts={"default": group_counts},
        sequence_map={"default": group_map},
        guide_gene_map={"default": guidemap},
    )
    try:
        model = chronos.Chronos(
            negative_control_sgrnas={"default": negative_controls},
            initial_screen_delay=initial_screen_delay,
            **common,
        )
    except Exception as first_error:
        logger.info(
            f"Chronos rejected the negative controls for {control_sample} "
            f"({first_error}); retrying without them."
        )
        try:
            model = chronos.Chronos(**common)
        except Exception as second_error:
            logger.warning(
                f"Chronos failed for {control_sample}: {second_error}. Skipping."
            )
            return

    model.train()
    os.makedirs(group_prefix, exist_ok=True)
    model.save(group_prefix, overwrite=True)


#### Analysis registry ####


def dry_function(*args, **kwargs):
    """Stand-in used by --dry-run so that the options are validated but nothing runs."""
    return None


ANALYSIS_FUNCTIONS = {
    "mageck": call_mageck_batch,
    "drugz": call_drugz_batch,
    "chronos": call_chronos_batch,
    "dry": dry_function,
}

#: The methods a configuration file may name.
AVAILABLE_ANALYSES = ("mageck", "drugz", "chronos")


#### Configuration ####


def process_control_map(controls: dict, samples) -> dict:
    """Expand the ALL and EXCEPT keywords into explicit sample lists."""
    controls = copy(controls)
    for control_map in controls.values():
        for control, treat_samples in control_map.items():
            if isinstance(treat_samples, str):
                treat_samples = [treat_samples]
            if treat_samples == ["ALL"]:
                treat_samples = [s for s in samples if s != control]
            elif isinstance(treat_samples, dict):
                treat_samples = [
                    s for s in samples if s not in treat_samples["EXCEPT"]
                ]
            control_map[control] = treat_samples
    return controls


def validate_required_arguments(arguments: dict) -> None:
    """Raise if the configuration cannot be run. Checks everything before
    reporting, so all problems surface at once."""
    problems: List[str] = []

    missing = [
        option
        for option in ("sample_reps", "experiment_id", "analyses", "file_prefix")
        if not arguments.get(option)
    ]
    if missing:
        problems.append(f"Missing or empty required option(s): {', '.join(missing)}")

    analyses = arguments.get("analyses") or []
    control_groups = arguments.get("control_groups") or {}

    for analysis in analyses:
        if "method" not in analysis:
            problems.append(f"Analysis has no method: {analysis}")

    expected_types = {
        "method": str, "kwargs": dict, "groups": list, "counts_file": str, "label": str,
    }
    for analysis in analyses:
        for option, expected in expected_types.items():
            if option in analysis and not isinstance(analysis[option], expected):
                problems.append(
                    f"Analysis option {option} should be {expected.__name__}, got "
                    f"{type(analysis[option]).__name__}: {analysis}"
                )

    samples = list(arguments.get("sample_reps") or {})
    hyphenated = [s for s in samples if "-" in s]
    if hyphenated:
        problems.append(
            f"Sample name(s) contain '-', which separates the two halves of a "
            f"comparison in output filenames: {', '.join(hyphenated)}"
        )

    unknown_samples = set()
    for control_map in control_groups.values():
        for control, treat_samples in control_map.items():
            if control not in samples:
                unknown_samples.add(control)
            for sample in treat_samples:
                if sample not in samples:
                    unknown_samples.add(sample)
    if unknown_samples:
        problems.append(
            "Sample(s) named in Control groups but absent from Sample details: "
            + ", ".join(sorted(unknown_samples))
        )

    unknown_groups = {
        group
        for analysis in analyses
        for group in analysis.get("groups", [])
        if group not in control_groups
    }
    if unknown_groups:
        problems.append(
            f"Group(s) named in Analyses but absent from Control groups: "
            f"{', '.join(sorted(unknown_groups))}. Control groups are: "
            f"{', '.join(control_groups)}"
        )

    if problems:
        for problem in problems:
            logger.error(problem)
        raise PipelineOptionsError(
            f"{len(problems)} problem(s) in the configuration. See above."
        )


def _check_counts_files(arguments: dict) -> None:
    """Confirm each counts file is tab separated and holds every replicate."""
    replicates = [
        replicate
        for group in arguments["sample_reps"].values()
        for replicate in group
    ]

    counts_files = set()
    global_counts = arguments.get("counts_file")
    if global_counts and os.path.isfile(global_counts):
        counts_files.add(global_counts)
    else:
        for analysis in arguments["analyses"]:
            if "counts_file" not in analysis:
                raise PipelineOptionsError(
                    f"No counts file for analysis {analysis.get('name', analysis)}. "
                    "Set one on the Analyses sheet or pass --counts."
                )
            counts_files.add(analysis["counts_file"])
            if analysis["method"] not in AVAILABLE_ANALYSES:
                raise PipelineOptionsError(
                    f"Unknown method {analysis['method']!r}. Available methods: "
                    f"{', '.join(AVAILABLE_ANALYSES)}."
                )

    if not counts_files:
        raise PipelineOptionsError("No counts file set for any analysis.")

    for fn in counts_files:
        fn = maybe_its_gz(fn)
        opener = gzip.open if fn.endswith(".gz") else open
        with opener(fn, "rt") as f:
            header = next(f)

        if "\t" not in header:
            raise PipelineOptionsError(
                f"No tabs in the first line of {fn}. Is it comma separated?\n\t{header}"
            )

        columns = header.strip().split("\t")
        missing = [r for r in replicates if r not in columns]
        if missing:
            raise PipelineOptionsError(
                f"Replicate(s) named in the workbook but absent from {fn}: "
                f"{', '.join(missing)}\nCounts columns: {', '.join(columns)}"
            )


def _filter_by_name(arguments: dict, key: str, target: str, describe: str) -> None:
    """Apply --run-groups / --run-analyses, which prune what gets run."""
    wanted = arguments.pop(key, None)
    if wanted is None:
        return
    if isinstance(wanted, str):
        wanted = [w.strip() for w in wanted.split(",")]

    if target == "control_groups":
        removed = [g for g in arguments["control_groups"] if g not in wanted]
        for group in removed:
            del arguments["control_groups"][group]
        remaining = list(arguments["control_groups"])
    else:
        kept = [a for a in arguments["analyses"] if a.get("name") in wanted]
        removed = [
            a.get("name") for a in arguments["analyses"] if a.get("name") not in wanted
        ]
        arguments["analyses"] = kept
        remaining = [a.get("name") for a in kept]

    if removed:
        logger.info(f"Removed {describe}: {', '.join(sorted(map(str, removed)))}")
    if not remaining:
        raise PipelineOptionsError(
            f"No {describe} left after applying --{key.replace('_', '-')} "
            f"{','.join(wanted)}."
        )
    logger.info(f"Running with {describe}: {', '.join(map(str, remaining))}")


def process_arguments(arguments: dict, delete_unrequired_args=True) -> dict:
    """Validate the configuration and expand its shorthand into what
    run_analyses expects."""
    validate_required_arguments(arguments)

    arguments["control_groups"] = process_control_map(
        arguments["control_groups"], list(arguments["sample_reps"])
    )

    # Results are filed under the experiment and, if set, the analysis version.
    arguments["output_dir"] = os.path.join(
        arguments["output_dir"],
        arguments["experiment_id"],
        arguments.get("analysis_version") or "",
    )
    logger.info(f"Output directory: {arguments['output_dir']}")

    for sample, replicates in arguments["sample_reps"].items():
        arguments["sample_reps"][sample] = list_not_str(replicates)

    _check_counts_files(arguments)
    _filter_by_name(arguments, "run_groups", "control_groups", "control groups")
    _filter_by_name(arguments, "run_analyses", "analyses", "analyses")

    if delete_unrequired_args:
        for key in ("config_file", "notes", "experiment_id", "analysis_version",
                    "labels"):
            arguments.pop(key, None)

    return arguments


def load_configuration_file(config_filename: str, counts_dir=".") -> dict:
    """Read an .xlsx workbook or a .json configuration."""
    if config_filename.endswith("xlsx"):
        return AnalysisWorkbook(config_filename, counts_dir=counts_dir).expd

    if config_filename.endswith("json"):
        # Lines beginning with # are treated as comments, which JSON has no
        # syntax for.
        with open(config_filename) as f:
            body = "\n".join(
                line for line in f if line.replace(" ", "")[:1] != "#"
            )

        def reject_duplicates(pairs):
            seen = {}
            for key, value in pairs:
                if key in seen:
                    raise ValueError(f"Duplicate key in JSON: {key!r}")
                seen[key] = value
            return seen

        return json.loads(body, object_pairs_hook=reject_duplicates)

    raise PipelineOptionsError(
        f"Unrecognised configuration file type: {config_filename}. "
        "Only .xlsx and .json are accepted."
    )


#### Running ####


def run_analyses(output_dir, file_prefix, sample_reps: Dict[str, list],
                 control_groups: Dict[str, Dict[str, list]], analyses: List[dict],
                 days_grown: Optional[Dict[str, list]] = None,
                 cell_line_hash: Optional[Dict[str, list]] = None,
                 counts_file=None, methods_kwargs: Optional[Dict] = None,
                 dont_log=False, overwrite=False, compjoiner=FILENAME_JOINER,
                 notes="", skip_method=None, dry_run=False, **unexpected) -> None:
    """Run every analysis named in `analyses`, then tabulate the results.

    days_grown and cell_line_hash default to empty because only Chronos needs
    them; a JSON configuration that omits them used to raise a TypeError here.
    """
    if unexpected:
        logger.warning(f"Ignoring unexpected arguments: {sorted(unexpected)}")
    if notes:
        logger.info(f"Notes: {notes}")

    days_grown = days_grown or {}
    cell_line_hash = cell_line_hash or {}

    if isinstance(skip_method, str):
        skip_method = [m.strip() for m in skip_method.split(",")]
    skip_method = list(skip_method or [])

    output_dir = str(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    if not dont_log:
        stamp = datetime.datetime.now().strftime("%Y-%m-%d_%Hh%Mm%Ss")
        logger.add(str(Path(output_dir, f"{file_prefix}log_{stamp}.txt")), level="INFO")

    logger.info(f"Full output directory: {os.path.realpath(output_dir)}")

    methods_used = {analysis["method"] for analysis in analyses}
    os.makedirs(Path(output_dir, "tables"), exist_ok=True)
    for method in methods_used:
        os.makedirs(Path(output_dir, method, "files"), exist_ok=True)

    # (method, results_prefix, table_prefix) -> the control map it was run with,
    # collected here so that tabulation happens once, after everything has run.
    # It used to sit inside the analysis loop and reuse whichever control map had
    # last leaked out of the inner loop.
    ran: Dict[tuple, dict] = {}

    for analysis in analyses:
        method = "dry" if dry_run else analysis["method"]

        if method in skip_method:
            logger.info(f"Skipping method {method}")
            continue

        # Defaults for a method can be set at the top level, and overridden per
        # analysis.
        default_kwargs = (methods_kwargs or {}).get(method) or {}
        kwargs = analysis.get("kwargs") or default_kwargs
        pseudocount = analysis.get("pseudocount", 1)
        groups = list_not_str(analysis.get("groups") or list(control_groups))

        for group in groups:
            control_map = control_groups[group]
            current_counts = maybe_its_gz(counts_file or analysis["counts_file"])
            out_prefix = str(Path(output_dir, method, "files", file_prefix))

            logger.info(
                f"Running {method} on group {group} with kwargs {kwargs} "
                f"and counts {current_counts}"
            )
            ran[(method, out_prefix, file_prefix)] = control_map

            ANALYSIS_FUNCTIONS[method](
                sample_reps, days_grown, cell_line_hash, control_map,
                current_counts, out_prefix, kwargs, pseudocount=pseudocount,
                overwrite=overwrite,
            )

    if dry_run:
        logger.info("Dry run: options validated, no analysis was run.")
        return

    for (method, results_prefix, table_prefix), control_map in ran.items():
        table = TABULATORS[method](results_prefix, cell_line_hash, control_map,
                                   compjoiner)
        table_fn = os.path.join(
            output_dir, "tables", f"{table_prefix}.{method}_table.csv"
        )
        logger.info(f"Writing table: {table_fn}")
        table.to_csv(table_fn, encoding="utf-8-sig")
