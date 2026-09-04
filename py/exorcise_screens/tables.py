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

from exorcise_screens.util import ARROW, iter_files_by_prefix, maybe_its_gz

MAGECK_STATISTICS = [
    "lfc", "fdr", "fdr_log10", "p", "p_log10", "pos_p", "neg_p", "neg_fdr", "pos_fdr",
]


def _empty_table() -> pd.DataFrame:
    return pd.DataFrame()


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

    columns = pd.MultiIndex.from_product([sorted(tables), tab.columns])
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

        tab = read_chronos(os.path.join(parent, fn, "gene_effect.hdf5"))
        tab = tab.reset_index()
        tab["index"] = tab["index"].astype(str)
        tab = tab.merge(hash_table, how="left").drop("index", axis=1)
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
            table[column] = sub_table[column]

    return table


#: Which tabulator belongs to which analysis method.
TABULATORS = {
    "mageck": tabulate_mageck,
    "drugz": tabulate_drugz,
    "chronos": tabulate_chronos,
}
