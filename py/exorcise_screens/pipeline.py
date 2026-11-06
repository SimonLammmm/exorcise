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

#: Guides matching these are treated as negative controls by Chronos: guides
#: expected to produce no phenotype, which is what it estimates excess variance
#: from and what it normalises readcounts against. The first pattern excludes
#: Exorcise's own exo_Non-targeting labels so that a library designer's declared
#: controls take precedence; the second accepts them.
#:
#: Both carried `.*_OR\d.*` until 3.1.6, inherited from crispr_tools, for
#: libraries that name their cutting controls after olfactory receptors.
#: Exorcise writes the gene symbol into exo_id, so against an exorcised library
#: that pattern matched every guide targeting a real olfactory receptor gene --
#: OR4F3, OR51B5, OR2I1P and the rest -- and handed Chronos thousands of
#: targeting guides as "expected no phenotype". No library in use here declares
#: OR cutting controls, so the pattern only ever fired as a false positive.
#: Tried in this order, and the first set Chronos actually accepts is the one
#: used. They are kept apart rather than combined because they are different
#: claims about a guide, and Chronos does not accept them interchangeably: a
#: declared non-targeting guide binds nothing, a cutting control cuts somewhere
#: harmless, and an Exorcise label means only that the guide did not map. A
#: screen whose cutting controls Chronos rejects may still fit on its
#: non-targeting ones, which is why a rejection moves to the next set rather
#: than straight to no controls at all.
CHRONOS_CONTROL_PATTERNS = (
    (r".*(?<!exo_)Non-targeting.*", "non-targeting controls declared by the library design"),
    (r".*_Cutting\d.*", "cutting controls declared by the library design"),
    (r".*Non-targeting.*", "Exorcise's unmapped-guide labels"),
)

#: Chronos raises rather than fit when fewer than 20 negative control sgRNAs
#: survive its own filtering, so a set smaller than this is worse than none: it
#: costs a fit that would otherwise have run on Chronos's default excess
#: variance. Such a set is passed over in favour of the next pattern.
CHRONOS_MIN_NEGATIVE_CONTROLS = 20


class PipelineOptionsError(Exception):
    """A configuration file problem that stops the pipeline before it starts."""


def read_counts(counts_file: str, index_col=None) -> pd.DataFrame:
    """Read a counts table, keeping the guide and gene columns as text.

    The first two columns of a counts file identify the guide and the gene. Most
    libraries name them something obviously textual, but a few number their
    guides 1, 2, 3... and at least one numbers its genes too. Pandas types a
    column of bare numbers as int64, and integer identifiers break things in
    three separate places:

    * Chronos calls ``.encode("utf8")`` on the column labels of the readcounts
      frame when it writes its HDF5 files, which raises AttributeError on an
      int after the model has already been fitted;
    * the negative controls are matched as strings, so they never line up with
      an integer guidemap and Chronos loses its excess-variance estimate;
    * the pseudocount is added to every non-object column, which means it lands
      on the identifiers as well as on the counts.

    Forcing the two identifier columns to str at the point of reading fixes all
    three. It is done by name rather than by position so that it applies whether
    or not the caller asks for the guide column as the index.
    """
    path = maybe_its_gz(counts_file)
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    as_text = {name: str for name in header[:2]}
    return pd.read_csv(path, sep="\t", index_col=index_col, dtype=as_text)


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
        counts = read_counts(counts_file, index_col=0)

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
    counts = read_counts(counts_file, index_col=0)
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
    counts = read_counts(counts_file)
    guide_column, gene_column = counts.columns[0], counts.columns[1]

    # Candidates are chosen from the labels, which a pseudocount does not touch,
    # but how well they were sequenced has to be read before one is added: a
    # pseudocount turns an undetected guide into a guide that reads exactly the
    # pseudocount, which is indistinguishable from a real low count by the time
    # Chronos sees it.
    control_candidates = _chronos_negative_control_candidates(
        counts[guide_column], counts[gene_column]
    )
    _report_control_detection(counts, guide_column, control_candidates, counts_file)

    if pseudocount > 1:
        numeric = counts.dtypes != object
        counts.loc[:, numeric] += pseudocount

    readcounts = (
        counts.drop(gene_column, axis=1)
        .set_index(guide_column)
        .transpose()
        .reset_index()
        .rename(columns={"index": "sequence_ID"})
    )

    guidemap = pd.DataFrame({"sgrna": counts[guide_column], "gene": counts[gene_column]})
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
            group_map, readcounts, guidemap, control_candidates, group_prefix,
            min(timepoints), control_sample,
        )


def _chronos_negative_control_candidates(guides: pd.Series, genes: pd.Series):
    """Ordered candidate sets of negative control guides, best claim first.

    Chronos uses these to estimate how much of the variance in a screen is noise
    rather than biology, and to normalise readcounts against. Without them it
    falls back to a flat default, which is a worse model, so it is worth looking
    properly - and worth having more than one answer ready, because Chronos
    checks the set it is given and rejects one it cannot normalise with.

    Which column names a control depends on the library. Exorcise writes the
    symbol into the guide ID, so `exo_id` carries it; libraries that were counted
    against their own design usually have opaque guide IDs like "ID_1" and name
    the control in the gene column instead. Both are searched, because only
    looking at the guide ID silently found nothing for the latter.

    Sets smaller than CHRONOS_MIN_NEGATIVE_CONTROLS are not offered at all:
    Chronos refuses to fit with them, so offering one only costs a fit that
    would have run on its default excess variance.
    """
    guides = guides.astype(str)
    genes = genes.astype(str)

    candidates = []
    seen = set()
    undersized = []

    for pattern, description in CHRONOS_CONTROL_PATTERNS:
        selected = guides[guides.str.match(pattern) | genes.str.match(pattern)]
        if len(selected) < CHRONOS_MIN_NEGATIVE_CONTROLS:
            if len(selected):
                undersized.append(f"{len(selected)} {description}")
            continue
        # A later pattern often matches everything an earlier one did. Offering
        # the same set twice would just repeat a rejection.
        fingerprint = frozenset(selected)
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        candidates.append((description, selected))

    if candidates:
        logger.info(
            "Negative control guides available to Chronos: "
            + "; ".join(f"{len(s)} from {d}" for d, s in candidates)
            + ". The first set Chronos accepts will be used."
        )
    elif undersized:
        logger.warning(
            "Too few negative control guides to fit with (found "
            + "; ".join(undersized)
            + f", and Chronos needs at least {CHRONOS_MIN_NEGATIVE_CONTROLS}). "
            "Continuing without them, on Chronos's default excess variance."
        )
    else:
        logger.warning(
            "No negative control guides were found, so Chronos cannot estimate "
            "excess variance and will fall back to its default. Controls are "
            "recognised by matching any of "
            + ", ".join("/" + p + "/" for p, _ in CHRONOS_CONTROL_PATTERNS)
            + " against either the guide or the gene column of the counts file."
        )
    return candidates


def _report_control_detection(counts, guide_column, candidates, counts_file) -> None:
    """Say how much of each candidate control set the sequencing actually saw.

    Chronos normalises each replicate so that the median of its negative
    controls lines up with the median in the pDNA samples. A control guide that
    was never sequenced carries nothing for that median, and a set that is
    mostly undetected gives a median that is really just the detection floor --
    or, once a pseudocount has been added, the pseudocount itself. Chronos then
    rejects the set, and what it says when it does is a paragraph about
    normalisation that never mentions the counts.

    So the state of the counts is worth reporting here, before any of that, in
    terms of the file the user can go and look at.
    """
    numeric = counts.select_dtypes("number")
    if numeric.empty:
        return

    guides = counts[guide_column].astype(str)
    for description, control_guides in candidates:
        rows = guides.isin(set(control_guides))
        block = numeric.loc[rows]
        if block.empty:
            continue
        detected = int((block > 0).all(axis=1).sum())
        if detected >= CHRONOS_MIN_NEGATIVE_CONTROLS:
            continue
        logger.warning(
            f"Only {detected} of the {len(block)} negative control guide(s) "
            f"from {description} were detected in every sample of "
            f"{os.path.basename(counts_file)}. Chronos normalises against the "
            "median of these, so a set this sparse carries almost no signal and "
            "it will probably reject the set. This is a property of the "
            "sequencing depth, not of the annotation."
        )


def _first_present(values):
    """The first usable value in a per-replicate list, or None.

    Days grown and Trajectory describe the sample, so every replicate of a
    sample carries the same value. Some workbooks fill them in on only one of
    the replicate rows, which is why this takes the first rather than requiring
    all of them.
    """
    for v in values or ():
        if v is None:
            continue
        text = str(v).strip()
        if text and text.lower() not in ("nan", "none"):
            return text
    return None


def _chronos_sequence_map(sample_reps, days_grown, cell_line_hash) -> pd.DataFrame:
    """The replicate/timepoint/trajectory table Chronos wants, one row per
    replicate.

    Chronos models counts at the replicate level, so every replicate has to be
    listed. Until 3.1.4 this was built by merging three per-replicate frames on
    `sample` and then calling `drop_duplicates(subset=["sample"])`, which kept
    exactly one replicate per sample and silently discarded the rest: a screen
    with three replicates handed Chronos a third of its data. The merge was also
    a cross product, so the surviving row's day and trajectory came from an
    arbitrary pairing rather than from that replicate.

    Rows whose day cannot be read are dropped, as before, because Chronos needs
    a timepoint for every sequence it is given.
    """
    rows = []
    dropped = 0

    for sample, replicates in sample_reps.items():
        day_text = _first_present(days_grown.get(sample))
        line = _first_present(cell_line_hash.get(sample))

        try:
            day = float(day_text)
        except (TypeError, ValueError):
            day = float("nan")

        for replicate in list_not_str(replicates) or ():
            if replicate is None or not str(replicate).strip():
                continue
            if math.isnan(day):
                dropped += 1
                continue
            rows.append({
                "sequence_ID": str(replicate),
                "sample": sample,
                "days": day,
                "cell_line_name": str(line),
                "pDNA_batch": 0,
            })

    if dropped:
        logger.warning(
            f"{dropped} replicate(s) have no readable Days grown and were left "
            "out of the Chronos sequence map."
        )

    sequence_map = pd.DataFrame(
        rows, columns=["sequence_ID", "sample", "days", "cell_line_name", "pDNA_batch"]
    )
    logger.info(
        f"Chronos sequence map: {len(sequence_map)} replicate(s) across "
        f"{sequence_map['sample'].nunique()} sample(s)."
    )
    return sequence_map


#: What Chronos itself uses for excess variance when it has no negative controls
#: to estimate it from. Supplying it explicitly is a workaround: some versions
#: leave it as an empty dict on this path and then reject their own value with
#: "excess_variance was passed as dict without key for 'default'".
CHRONOS_DEFAULT_EXCESS_VARIANCE = 0.02


def _fit_chronos_without_controls(chronos, common, initial_screen_delay,
                                  control_sample):
    """Fit a model for a screen with no usable negative controls.

    Returns None if it cannot be fitted, having said why.
    """
    # Excess variance has to be supplied explicitly on this path. Chronos
    # estimates it from negative control sgRNAs, and there are none here, so it
    # leaves the value an empty dict and then rejects it with "excess_variance
    # was passed as dict without key for 'default'". Until 3.1.6 two of the
    # three attempts here omitted it and so could never have succeeded; they
    # only added that message to every failure report.
    attempts = (
        ("supplying the default excess variance",
         dict(excess_variance=CHRONOS_DEFAULT_EXCESS_VARIANCE,
              initial_screen_delay=initial_screen_delay)),
        # A Chronos too old to accept initial_screen_delay. The excess variance
        # is not optional, so that is what the retry drops.
        ("supplying the default excess variance, without initial_screen_delay",
         dict(excess_variance=CHRONOS_DEFAULT_EXCESS_VARIANCE)),
    )

    failures = []
    for description, extra in attempts:
        try:
            return chronos.Chronos(**common, **extra)
        except TypeError as error:
            # An argument this version of Chronos does not accept. Worth trying
            # the next, simpler, form.
            failures.append(f"{description}: {error}")
        except Exception as error:
            failures.append(f"{description}: {error}")

    # Attribute the failure to the right thing. "estimated efficacy outside
    # bounds" is Chronos rejecting the screen's own fold changes -- it found no
    # depletion to scale against in at least one replicate -- and no amount of
    # relabelling negative controls will change that.
    rejected_data = any("efficacy outside bounds" in f for f in failures)
    if rejected_data:
        advice = (
            "\n  Chronos rejected the screen's fold changes rather than its "
            "configuration: in at least one replicate the guides that deplete "
            "most across the group do not deplete at all, so it cannot estimate "
            "how efficiently that replicate was edited. That is usually a "
            "vehicle or untreated arm carried as its own trajectory, a screen "
            "too short for depletion to show, or a control timepoint whose "
            "replicates disagree. MAGeCK and DrugZ are unaffected."
        )
    else:
        advice = (
            "\n  Negative controls are what Chronos estimates excess variance "
            "from; check that the counts file labels them in either its guide "
            "or its gene column."
        )

    logger.warning(
        f"Chronos could not fit {control_sample} without negative controls, so "
        f"this comparison has no Chronos result. Attempts:\n    "
        + "\n    ".join(failures)
        + advice
    )
    return None


def _mask_outgrowths(chronos, group_counts, group_map, guidemap,
                     control_sample) -> None:
    """Blank out clonal outgrowths before fitting, when that is safe to do.

    Chronos's nan_outgrowths looks for guides whose late counts run away from
    the rest, comparing each guide against the other guides for the same gene.
    To do that it needs something to compare against, and it divides by a count
    of comparable cases. In a one-guide-per-gene library where every trajectory
    also has a single late replicate, there is no such case anywhere in the
    screen and it raises ZeroDivisionError.

    A library with several guides per gene always survives that filter however
    few replicates it has. A one-guide-per-gene library does not, so both
    conditions have to hold before masking is skipped - an earlier version of
    this guard keyed on the replicate count alone and would have switched
    masking off for the great majority of screens, which do not need it off.
    """
    late = group_map.loc[group_map["cell_line_name"] != "pDNA"]
    per_trajectory = late.groupby("cell_line_name").size()

    try:
        per_gene = guidemap.groupby("gene").size()
        one_guide_per_gene = len(per_gene) > 0 and int(per_gene.max()) <= 1
    except Exception:
        one_guide_per_gene = False

    if one_guide_per_gene and len(per_trajectory) and (per_trajectory <= 1).all():
        logger.info(
            f"Not masking outgrowths for {control_sample}: the library has one "
            "guide per gene and every trajectory has a single late replicate, "
            "so there is nothing to compare a guide against."
        )
        return

    try:
        chronos.nan_outgrowths(group_counts, group_map, guidemap)
    except ZeroDivisionError:
        logger.warning(
            f"Chronos found no comparable log fold changes while masking "
            f"outgrowths for {control_sample}, so none were masked. The fit "
            "continues on the unmasked counts."
        )


def _fit_chronos(group_map, readcounts, guidemap, control_candidates, group_prefix,
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
    _mask_outgrowths(chronos, group_counts, group_map, guidemap, control_sample)

    common = dict(
        readcounts={"default": group_counts},
        sequence_map={"default": group_map},
        guide_gene_map={"default": guidemap},
    )

    # Each candidate set in turn. Chronos checks the controls it is given --
    # it has to normalise the readcounts so their medians line up against the
    # pDNA samples -- and a set it rejects says nothing about the next one, so a
    # rejection moves on rather than abandoning controls altogether.
    model = None
    for description, negative_controls in control_candidates:
        try:
            model = chronos.Chronos(
                negative_control_sgrnas={"default": negative_controls},
                initial_screen_delay=initial_screen_delay,
                **common,
            )
            logger.info(
                f"Chronos fitted {control_sample} on {len(negative_controls)} "
                f"negative control guide(s) from {description}."
            )
            break
        except Exception as error:
            logger.warning(
                f"Chronos rejected the {len(negative_controls)} negative "
                f"control guide(s) from {description} for {control_sample} "
                f"({error})."
            )

    if model is None:
        model = _fit_chronos_without_controls(
            chronos, common, initial_screen_delay, control_sample
        )
    if model is None:
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

#: Methods this pipeline computes, each from a counts table.
COUNTS_ANALYSES = ("mageck", "drugz", "chronos")

#: Methods whose gene-level results are produced outside this pipeline and
#: supplied as a finished table. Nothing is computed for them and no counts
#: table is needed; the pipeline validates what was supplied and leaves it
#: alone. A workbook using one of these exists to carry the metadata that the
#: database build needs.
SUPPLIED_ANALYSES = ("manual",)

#: Statistics every supplied table must carry for each comparison. `score` and
#: `pval` are what database.build._select_stat_columns reads back out.
SUPPLIED_STATISTICS = ("score", "pval")

#: The methods a configuration file may name.
AVAILABLE_ANALYSES = COUNTS_ANALYSES + SUPPLIED_ANALYSES


#### Supplied results ####


def expected_comparisons(control_groups: Dict[str, Dict[str, list]],
                         groups, compjoiner=FILENAME_JOINER) -> List[str]:
    """The comparison labels a set of control groups should produce.

    Built the same way the analysis functions name their outputs, including
    skipping a sample compared against itself, so that a supplied table can be
    held to the same naming the computed ones use.
    """
    labels: List[str] = []
    for group in groups:
        for control, treat_samples in control_groups[group].items():
            for treat in list_not_str(treat_samples):
                if treat == control:
                    continue
                labels.append(f"{control}{compjoiner}{treat}")
    return labels


def supplied_table_path(output_dir, file_prefix: str, method: str) -> str:
    """Where a supplied results table is expected to sit."""
    return maybe_its_gz(os.path.join(
        str(output_dir), "tables", f"{file_prefix}.{method}_table.csv"
    ))


def check_supplied_table(method: str, output_dir, file_prefix: str,
                         control_groups: Dict[str, Dict[str, list]], groups,
                         compjoiner=FILENAME_JOINER) -> str:
    """Validate a results table that was produced outside this pipeline.

    A supplied table is not checkable the way a computed one is: nothing here
    made it, so the only thing standing between a typo and a silently wrong
    database is this function. It is strict on purpose, and it reports every
    problem it finds rather than stopping at the first.

    The table is the double-headered CSV the tabulators write: comparison on
    the first header row, statistic on the second, genes down the index.
    """
    path = supplied_table_path(output_dir, file_prefix, method)
    if not os.path.isfile(path):
        raise PipelineOptionsError(
            f"The {method!r} method supplies its own results, but no table was "
            f"found at {path}. Nothing computes it, so it has to be put there "
            "before the pipeline runs."
        )

    try:
        table = pd.read_csv(path, index_col=0, header=[0, 1])
    except Exception as error:
        raise PipelineOptionsError(
            f"Could not read the supplied {method!r} table {path}: {error}. It "
            "should be a CSV with two header rows, comparison then statistic, "
            "and the gene in the first column."
        ) from error

    problems: List[str] = []

    # Unnamed: N_level_1 is what pandas calls a blank cell on the second header
    # row, which means the column has a comparison but no statistic.
    unlabelled = [
        str(c[0]) for c in table.columns if str(c[1]).startswith("Unnamed:")
    ]
    if unlabelled:
        problems.append(
            f"{len(unlabelled)} column(s) have no statistic on the second "
            f"header row, e.g. {', '.join(sorted(set(unlabelled))[:3])}"
        )

    found = [str(c[0]) for c in table.columns]
    expected = expected_comparisons(control_groups, groups, compjoiner)

    repeated = sorted({c for c in expected if expected.count(c) > 1})
    if repeated:
        problems.append(
            f"{len(repeated)} comparison(s) are defined more than once in "
            f"Control groups: {', '.join(repeated[:5])}"
        )

    missing = sorted(set(expected) - set(found))
    if missing:
        problems.append(
            f"{len(missing)} comparison(s) in the workbook are absent from the "
            f"table: {', '.join(missing[:5])}"
            + (" ..." if len(missing) > 5 else "")
        )

    extra = sorted(set(found) - set(expected))
    if extra:
        problems.append(
            f"{len(extra)} comparison(s) in the table are absent from the "
            f"workbook: {', '.join(extra[:5])}"
            + (" ..." if len(extra) > 5 else "")
        )

    # Every comparison needs the full set of statistics, and needs each of them
    # exactly once. A duplicate would be silently dropped on the way into the
    # database, taking an arbitrary one of the two.
    statistics: Dict[str, List[str]] = {}
    for comparison, statistic in table.columns:
        statistics.setdefault(str(comparison), []).append(str(statistic))

    incomplete, duplicated = [], []
    for comparison in set(found) & set(expected):
        present = statistics[comparison]
        if set(SUPPLIED_STATISTICS) - set(present):
            incomplete.append(
                f"{comparison} (has {', '.join(sorted(set(present)))})"
            )
        if len(present) != len(set(present)):
            duplicated.append(comparison)
    if incomplete:
        problems.append(
            f"{len(incomplete)} comparison(s) are missing one of "
            f"{', '.join(SUPPLIED_STATISTICS)}: {'; '.join(sorted(incomplete)[:3])}"
        )
    if duplicated:
        problems.append(
            f"{len(duplicated)} comparison(s) have a repeated statistic "
            f"column: {', '.join(sorted(duplicated)[:5])}"
        )

    if table.index.isna().any():
        problems.append(f"{int(table.index.isna().sum())} row(s) have no gene")
    if not table.index.is_unique:
        repeats = table.index[table.index.duplicated()].unique()
        problems.append(
            f"{len(repeats)} gene(s) appear on more than one row: "
            f"{', '.join(map(str, repeats[:5]))}"
        )
    if not len(table.index):
        problems.append("the table has no genes")

    if problems:
        for problem in problems:
            logger.error(problem)
        raise PipelineOptionsError(
            f"{len(problems)} problem(s) in the supplied {method!r} table "
            f"{path}. See above."
        )

    logger.info(
        f"Supplied {method!r} table {path} checked: {len(set(found))} "
        f"comparison(s) matching the workbook, {len(table.index)} gene(s)."
    )
    return path


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

    # A '-' in a sample name is ambiguous in the comparison labels that the
    # analysis functions build their output filenames from, and the tabulators
    # rewrite the joiner with a plain str.replace, which would eat a hyphen
    # inside a name. Neither of those happens for a method that supplies its
    # own results: its table is matched by exact string and is never rewritten.
    # So the rule is applied to the samples that a computed analysis could
    # actually name, rather than to every sample in the workbook.
    computed_groups = {
        group
        for analysis in analyses
        if analysis.get("method") in COUNTS_ANALYSES
        for group in (analysis.get("groups") or list(control_groups))
    }
    computed_samples = {
        sample
        for group in computed_groups
        if group in control_groups
        for control, treat_samples in control_groups[group].items()
        for sample in [control, *list_not_str(treat_samples)]
    }
    hyphenated = sorted(s for s in computed_samples if "-" in s)
    if hyphenated:
        problems.append(
            f"{len(hyphenated)} sample name(s) contain '-', which separates "
            f"the two halves of a comparison in output filenames: "
            f"{', '.join(hyphenated[:5])}"
            + (f" ... and {len(hyphenated) - 5} more" if len(hyphenated) > 5 else "")
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


def _check_methods(arguments: dict) -> None:
    """Reject method names nothing knows how to handle.

    This used to live inside _check_counts_files, in the branch taken only when
    no counts file was given on the command line, so an unknown method passed
    validation whenever --counts named a file.
    """
    unknown = sorted({
        analysis["method"]
        for analysis in arguments["analyses"]
        if analysis["method"] not in AVAILABLE_ANALYSES
    })
    if unknown:
        raise PipelineOptionsError(
            f"Unknown method(s) {', '.join(repr(m) for m in unknown)}. "
            f"Available methods: {', '.join(AVAILABLE_ANALYSES)}."
        )


def _check_counts_files(arguments: dict) -> None:
    """Confirm each counts file is tab separated and holds every replicate.

    Only the methods that read counts are considered. A workbook whose
    analyses all supply their own results needs no counts table, and whatever
    its Counts file cell happens to say is ignored rather than resolved.
    """
    needs_counts = [
        analysis
        for analysis in arguments["analyses"]
        if analysis["method"] in COUNTS_ANALYSES
    ]
    if not needs_counts:
        supplied = sorted({a["method"] for a in arguments["analyses"]})
        logger.info(
            f"No counts table is required: every analysis supplies its own "
            f"results ({', '.join(supplied)})."
        )
        return

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
        for analysis in needs_counts:
            if "counts_file" not in analysis:
                raise PipelineOptionsError(
                    f"No counts file for analysis {analysis.get('name', analysis)} "
                    f"(method {analysis['method']}). Set one on the Analyses "
                    "sheet or pass --counts."
                )
            counts_files.add(analysis["counts_file"])

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

    _check_methods(arguments)
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
        # A supplied method writes nothing, so it gets no working directory.
        if method in SUPPLIED_ANALYSES:
            continue
        os.makedirs(Path(output_dir, method, "files"), exist_ok=True)

    # (method, results_prefix, table_prefix) -> the control map it was run with,
    # collected here so that tabulation happens once, after everything has run.
    # It used to sit inside the analysis loop and reuse whichever control map had
    # last leaked out of the inner loop.
    ran: Dict[tuple, dict] = {}

    for analysis in analyses:
        method = analysis["method"]

        if method in skip_method:
            logger.info(f"Skipping method {method}")
            continue

        groups = list_not_str(analysis.get("groups") or list(control_groups))

        # Nothing computes a supplied method, so the table that was put there
        # is checked instead and then left exactly as it is. Checking is
        # read-only, so it is worth doing on a dry run too - validating the
        # configuration is the whole point of one.
        if method in SUPPLIED_ANALYSES:
            check_supplied_table(
                method, output_dir, file_prefix, control_groups, groups,
                compjoiner,
            )
            continue

        if dry_run:
            method = "dry"

        # Defaults for a method can be set at the top level, and overridden per
        # analysis.
        default_kwargs = (methods_kwargs or {}).get(method) or {}
        kwargs = analysis.get("kwargs") or default_kwargs
        pseudocount = analysis.get("pseudocount", 1)

        for group in groups:
            control_map = control_groups[group]
            current_counts = maybe_its_gz(counts_file or analysis["counts_file"])
            out_prefix = str(Path(output_dir, method, "files", file_prefix))

            logger.info(
                f"Running {method} on group {group} with kwargs {kwargs} "
                f"and counts {current_counts}"
            )
            # Accumulated, not assigned. Several groups can run into one
            # prefix and tabulation happens once over the lot, so assigning
            # kept only the last group's map. That was harmless while nothing
            # read it; it is not now that tabulate_chronos restricts its
            # columns to what the map declares, and would have dropped every
            # earlier group's comparisons as undeclared.
            declared = ran.setdefault((method, out_prefix, file_prefix), {})
            for ctrl, treats in control_map.items():
                known = declared.setdefault(ctrl, [])
                for treat in list_not_str(treats):
                    if treat not in known:
                        known.append(treat)

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

        # An empty table means every comparison failed. Each failure was warned
        # about as it happened, but those scroll past in a long run and the file
        # that gets left behind looks like a result: the database build reads an
        # unparseable or empty table as "this method contributed nothing" rather
        # than as an error. Say so once, plainly, at the point it is known.
        if table.empty:
            logger.warning(
                f"{method} produced no results for any comparison, so "
                f"{table_fn} is empty. The reason for each comparison is in the "
                "warnings above. Nothing downstream treats this as an error, so "
                "this is the only place it is reported."
            )
