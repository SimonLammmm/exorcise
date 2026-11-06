"""Turn each analysis program's output files into one table per experiment.

Every tabulate_* function returns a DataFrame indexed by gene with a two-level
column index of (comparison, statistic), which is the shape the database builder
expects.

Was the tabulate_* half of crispr_tools/tools.py.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from loguru import logger

from exorcise_screens.util import (
    ARROW,
    iter_files_by_prefix,
    list_not_str,
    maybe_its_gz,
)

MAGECK_STATISTICS = [
    "lfc", "fdr", "fdr_log10", "p", "p_log10", "pos_p", "neg_p", "neg_fdr", "pos_fdr",
]


def _empty_table() -> pd.DataFrame:
    return pd.DataFrame()


def _attested_comparisons(ctrl_map, compjoiner, spell=None):
    """The comparison labels the control map attests, spelled as the tabulator
    spells the ones it finds on disk.

    Returns None when there is no map to check against, which is not the same as
    an empty set: None means nothing can be concluded, empty means the map
    attests nothing.

    `spell` exists because the tabulators do not name a comparison the same way.
    MAGeCK strips every dot out of the filename it parses, so a sample named
    `plate.1` reaches the table as `plate1`; DrugZ strips only the separator and
    keeps the rest. Comparing a declared label with one read off disk therefore
    has to go through whatever the reading tabulator did, or every comparison on
    a screen whose samples contain dots looks unattested.
    """
    if not ctrl_map:
        return None
    labels = set()
    for control, treats in ctrl_map.items():
        for treat in list_not_str(treats):
            if str(treat) == str(control):
                continue
            label = f"{control}{compjoiner}{treat}"
            labels.add(spell(label) if spell else label)
    return labels


def _drop_unattested(tables, ctrl_map, compjoiner, method, spell=None):
    """Keep only the comparisons the workbook's control map attests.

    Results accumulate on disk. Edit a workbook and re-run without clearing the
    results directory first and the files from the previous shape of the screen
    are still sitting there, named after comparisons the workbook no longer
    makes. The tabulators find their inputs by globbing the results prefix, so
    those files were tabulated too, and the stale comparisons reached the
    database as statistics with no metadata describing them - nothing names them
    in Experiment or Comparison details, because nothing in the workbook names
    them at all.

    So the control map decides. A comparison it does not attest is reported and
    left out rather than carried forward silently. Clearing the results
    directory before re-running remains the way to avoid the situation; this
    makes failing to do so visible instead of quiet.
    """
    attested = _attested_comparisons(ctrl_map, compjoiner, spell)
    if attested is None:
        return tables

    stale = sorted(c for c in tables if c not in attested)
    if not stale:
        return tables

    logger.warning(
        f"{method}: {len(stale)} comparison(s) found in the results directory "
        f"are not attested by the workbook's control map and have been left out "
        f"of the table: {', '.join(stale)}. These are results of a shape the "
        "workbook no longer describes, left behind by a re-run over an "
        "uncleared results directory. They have no metadata, so adding them "
        "would put statistics in the database that nothing describes. Clear the "
        "results directory and re-run to be rid of them."
    )
    return {c: tab for c, tab in tables.items() if c in attested}


def tabulate_mageck(prefix, cell_line_hash=None, ctrl_map=None,
                    compjoiner=ARROW) -> pd.DataFrame:
    """Collect MAGeCK gene_summary files into one table."""
    prefix = Path(prefix)
    tables = {}
    tab = None

    for fn in iter_files_by_prefix(prefix, ".gene_summary.txt"):
        mtab = pd.read_csv(prefix.parent / fn, sep="\t", index_col=0)
        tab = pd.DataFrame(index=mtab.index)

        # MAGeCK reports the sensitising and suppressing directions separately.
        tab.loc[:, "neg_p"] = mtab["neg|p-value"]
        tab.loc[:, "pos_p"] = mtab["pos|p-value"]
        tab.loc[:, "neg_fdr"] = mtab["neg|fdr"]
        tab.loc[:, "pos_fdr"] = mtab["pos|fdr"]
        tab.loc[:, "lfc"] = mtab.loc[:, "neg|lfc"]

        # Collapse the two directions into one, taking whichever applies to the
        # sign of the fold change.
        enriched = mtab["pos|lfc"] > 0
        for statistic in ("fdr", "p-value"):
            key = statistic.replace("-value", "")
            tab.loc[enriched, key] = mtab.loc[enriched, f"pos|{statistic}"]
            tab.loc[~enriched, key] = mtab.loc[~enriched, f"neg|{statistic}"]
            tab.loc[:, f"{key}_log10"] = -np.log10(tab[key])

        comparison = (
            fn.split(prefix.stem, 1)[1]
            .split(".gene_s")[0]
            .replace("-", compjoiner)
            .replace(".", "")
        )
        tables[comparison] = tab

    if tab is None:
        raise FileNotFoundError(
            f"No .gene_summary.txt files found with prefix {prefix}. Did MAGeCK run?"
        )

    # Dots are stripped from the parsed filename above, so the attested labels
    # have to lose theirs too before they can be compared.
    tables = _drop_unattested(tables, ctrl_map, compjoiner, "MAGeCK",
                              spell=lambda s: s.replace(".", ""))
    if not tables:
        return _empty_table()

    columns = pd.MultiIndex.from_product([sorted(tables), MAGECK_STATISTICS])
    table = pd.DataFrame(index=tab.index, columns=columns)

    for comparison, tab in tables.items():
        # A gene absent from one comparison gets the null result rather than NaN.
        absent = tab.isna().any(axis=1)
        tab.loc[absent, :] = 1
        tab.loc[absent, ["lfc", "fdr_log10", "p_log10"]] = 0
        table[comparison] = tab

    return table


def tabulate_drugz(prefix, cell_line_hash=None, ctrl_map=None,
                   compjoiner=ARROW) -> pd.DataFrame:
    """Collect DrugZ result files into one table."""
    prefix = Path(prefix)
    tables = {}
    tab = None

    for fn in iter_files_by_prefix(prefix):
        comparison = fn.split(prefix.stem, 1)[1]
        comparison = re.sub(r"\.tsv.*$", "", comparison)
        comparison = re.sub(r"^\.", "", comparison)
        comparison = comparison.replace("-", compjoiner)

        tab = pd.read_csv(os.path.join(prefix.parent, fn), sep="\t", index_col=0)
        tab.index.name = "gene"
        tab = tab.loc[:, ["normZ", "pval_synth", "fdr_synth", "pval_supp", "fdr_supp"]]
        tab.columns = ["normZ", "neg_p", "neg_fdr", "pos_p", "pos_fdr"]

        # Report whichever direction is more significant.
        for statistic in ("p", "fdr"):
            best = tab.loc[:, [f"neg_{statistic}", f"pos_{statistic}"]].min(axis=1)
            tab.loc[:, statistic] = best
            tab.loc[:, f"{statistic}_log10"] = -np.log10(best)

        tables[comparison] = tab

    if tab is None:
        raise FileNotFoundError(
            f"No DrugZ result files found with prefix {prefix}. Did DrugZ run?"
        )

    # DrugZ keeps the dots inside a sample name, so the labels compare directly.
    statistics = list(tab.columns)
    tables = _drop_unattested(tables, ctrl_map, compjoiner, "DrugZ")
    if not tables:
        return _empty_table()

    columns = pd.MultiIndex.from_product([sorted(tables), statistics])
    table = pd.DataFrame(index=tab.index, columns=columns)
    for comparison, tab in tables.items():
        table[comparison] = tab

    return table


def read_chronos(filename):
    """Read a Chronos gene_effect HDF5 model into a DataFrame.

    Chronos writes plain HDF5; if only a gzipped copy is present it has to be
    expanded first because h5py cannot read through gzip.
    """
    import h5py  # imported here so that MAGeCK-only runs need not have it

    gz = maybe_its_gz(filename)
    expanded = False
    if gz != filename:
        import gzip

        with gzip.open(gz, "rb") as g:
            content = g.read()
        with open(filename, "wb") as o:
            o.write(content)
        expanded = True

    try:
        with h5py.File(filename, "r") as f:
            data = list(f["data"])
            rows = [b.decode("utf-8") for b in f["dim_0"]]
            columns = [b.decode("utf-8") for b in f["dim_1"]]
    finally:
        if expanded:
            os.remove(filename)

    return pd.DataFrame(data).set_axis(rows, axis=0).set_axis(columns, axis=1)


def _declared_treatments(level, ctrl_map):
    """The samples this control was declared against, or None if not declared.

    None and the empty set mean different things: None is "the map does not
    cover this control, so nothing can be concluded", and an empty set is "the
    map covers it and names nothing".
    """
    if not ctrl_map:
        return None
    for control, treats in ctrl_map.items():
        if str(control) == str(level):
            return {str(t) for t in list_not_str(treats) if str(t) != str(level)}
    return None


def _restrict_to_declared(tab, level, ctrl_map):
    """Keep only the samples this control group actually declared.

    Chronos reports one row per trajectory, and a trajectory is one arm
    followed over time, so its Trajectory value is shared by every sample on
    it. Mapping those rows back to sample names therefore yields a column for
    every sample sharing a trajectory with one Chronos modelled, whatever the
    workbook asked for. That is how a control came to be compared against
    itself, how samples left out of every comparison appeared anyway, and how
    samples were paired with controls they were never declared against.

    `ctrl_map` is the comparison map the pipeline already assembled and has
    always passed in; until 3.1.9 every tabulator accepted it and none read it.
    Where it names this control, its treatment samples are the whole of what
    belongs in the table.

    Where it does not -- an older results directory, or a tabulation run
    without a map -- nothing is assumed beyond the one case that is wrong on
    its own terms: a control compared against itself.
    """
    samples = tab["sample"].astype(str)
    declared = _declared_treatments(level, ctrl_map)

    if declared is None:
        return tab.loc[samples != str(level)]

    keep = samples.isin(declared)
    dropped = sorted(set(samples[~keep]))
    if dropped:
        logger.info(
            f"Chronos: dropping {len(dropped)} undeclared column(s) from "
            f"{level}: {', '.join(dropped)}. These share a trajectory with a "
            "sample the model covers but were not declared against this "
            "control."
        )
    return tab.loc[keep]


def tabulate_chronos(prefix, cell_line_hash, ctrl_map=None,
                     compjoiner=ARROW) -> pd.DataFrame:
    """Collect Chronos gene effect models into one table.

    Chronos writes one directory per analysis, and identifies replicates by the
    trajectory hash rather than by sample name, so the hash has to be mapped back
    to sample names here.
    """
    prefix_path = Path(prefix)
    parent = str(prefix_path.parent) + os.sep
    stem = prefix_path.name + r"\."

    hash_table = (
        pd.DataFrame.from_dict(cell_line_hash, orient="index")
        .unstack()
        .reset_index()
        .drop("level_0", axis=1)
        .rename(columns={"level_1": "sample", 0: "index"})
        .drop_duplicates()
    )
    hash_table = hash_table[[i is not None for i in hash_table["index"]]]
    hash_table["index"] = hash_table["index"].astype(str)

    tables = {}
    tab = None
    for fn in iter_files_by_prefix(prefix_path):
        level = re.sub(stem + r"(.+)", r"\1", fn)

        gene_effect_path = os.path.join(parent, fn, "gene_effect.hdf5")
        # Through maybe_its_gz, because a finished result is routinely gzipped
        # afterwards and read_chronos expands one itself. Until 3.1.8 this
        # looked only for the plain name, so a gzipped model was skipped here
        # before read_chronos ever saw it: call_chronos_batch reported the
        # result as already present and skipped recomputing it, and the
        # tabulator then reported that no comparison had produced anything.
        # Gunzipping by hand made the same run work, which is the tell.
        if not os.path.exists(maybe_its_gz(gene_effect_path)):
            # Skip levels where Chronos failed to produce gene_effect.hdf5
            continue

        try:
            tab = read_chronos(gene_effect_path)
        except (FileNotFoundError, KeyError, OSError):
            # Skip levels where Chronos output is corrupted
            continue

        tab = tab.reset_index()
        tab["index"] = tab["index"].astype(str)
        tab = tab.merge(hash_table, how="left").drop("index", axis=1)
        tab = _restrict_to_declared(tab, level, ctrl_map)
        if tab.empty:
            continue
        tab = tab.set_index("sample")
        tab = tab.set_axis([level + compjoiner + str(c) for c in tab.index])
        tab = tab.transpose()
        tab.index.name = "gene"
        tables[level] = tab

    if not tables:
        return _empty_table()

    comparisons = sorted(
        column for table in tables.values() for column in table.keys()
    )
    columns = pd.MultiIndex.from_product([comparisons, ["chronos_score"]])
    table = pd.DataFrame(index=tab.index, columns=columns)
    for sub_table in tables.values():
        for column in sub_table.keys():
            # Reindex to match the master table's index to avoid duplicate label errors
            # when indices don't align perfectly across different Chronos results
            try:
                table.loc[:, column] = sub_table[column]
            except (ValueError, KeyError):
                # If assignment fails due to index mismatch, skip this column
                continue

    return table


#: Which tabulator belongs to which analysis method.
TABULATORS = {
    "mageck": tabulate_mageck,
    "drugz": tabulate_drugz,
    "chronos": tabulate_chronos,
}
