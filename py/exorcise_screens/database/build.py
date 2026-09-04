"""Build and update the screen results database.

Given the analysis workbooks for one or more experiments, plus the results and
counts directories the pipeline wrote, this produces a directory containing:

    database.db                     gene and per-comparison statistics
    experiments_metadata.csv.gz     one row per experiment
    comparisons_metadata.csv.gz     one row per comparison

Was crispr_screen_viewer/update_database.py. Nothing here touches the network,
and nothing imports Dash: the web viewer is not part of this package.
"""

from __future__ import annotations

import glob
import os
import sys
import typing
from dataclasses import dataclass
from pathlib import Path
from typing import Collection, Dict, List, Optional, Sequence, Type, Union

import numpy as np
import pandas as pd
import sqlalchemy
from sqlalchemy import Engine, create_engine
from sqlalchemy.orm import Session
from loguru import logger

from exorcise_screens.database.metadata import (
    ANALYSESTYPES,
    AnalysisType,
    DB_FILES,
    MetadataTables,
    get_db_url,
)
from exorcise_screens.database.schema import GeneTable, StatTable, TableBase
from exorcise_screens.util import (
    TIMEPOINT_LABELS,
    df_rename_columns,
    is_nt,
    is_temp_file,
    load_stats_csv,
    maybe_its_gz,
    normalise_text,
)
from exorcise_screens.workbook import AnalysisWorkbook

#: Separates a control from a treatment when a comparison is described in prose.
FAT_ARROW = "➤"

#: Labels that have been used for the same field across template revisions.
EXPERIMENT_COLUMN_RENAMER = {
    "Experiment name": "Experiment ID",
    "Analysis name": "Experiment ID",  # old label for experiment name
    "Experiment description (a few sentances)": "Experiment description",
    "Experiment description (a few sentences)": "Experiment description",
    # Reference information used to live in a field called Citation. Kept for
    # backwards compatibility, though it would clash if a real citation string
    # were ever put in Experiment details.
    "Citation": "Reference",
    "Date screen completed (yyyy-mm-dd)": "Date",
    "Date screen completed": "Date",
    "Date published": "Date",
}


class ConfirmationRequired(Exception):
    """A destructive action needed consent that could not be obtained."""


@dataclass
class AnalysisInfo:
    """One experiment: its workbook, and where its files are."""

    experiment_id: str
    analysis_workbook: AnalysisWorkbook
    counts_path: str
    results_paths: Dict[AnalysisType, Union[Path, str]]


#### Locating inputs ####


def get_paths(details_xlsx: List[str], results_dir: Union[str, Path],
              count_dir: Union[str, Path],
              analysis_filename_prefix: str = "result") -> List[AnalysisInfo]:
    """Resolve each workbook to its counts file and result tables.

    Assumes the flat layout the pipeline writes: results under
    <results_dir>/<experiment_id>/tables/.
    """
    infos: List[AnalysisInfo] = []

    for fn in details_xlsx:
        if is_temp_file(fn):
            logger.debug(f"Skipping temporary file {fn}")
            continue

        workbook = AnalysisWorkbook(fn)
        experiment_id = workbook.experiment_details["Experiment name"]

        counts_files = workbook.analyses["Counts file"].dropna().unique()
        counts_files = [c for c in counts_files if str(c).strip()]
        if len(counts_files) != 1:
            raise RuntimeError(
                f"{fn} names {len(counts_files)} counts files "
                f"({', '.join(map(str, counts_files))}). Exactly one is supported."
            )

        counts_path = maybe_its_gz(os.path.join(count_dir, counts_files[0]))
        if not os.path.isfile(counts_path):
            raise FileNotFoundError(
                f"Counts file for experiment {experiment_id} not found at "
                f"{counts_path}. Is --counts-dir right?"
            )

        results_paths = {}
        for analysis_type in ANALYSESTYPES:
            table_fn = maybe_its_gz(os.path.join(
                results_dir, experiment_id, "tables",
                f"{analysis_filename_prefix}.{analysis_type.name}_table.csv",
            ))
            if os.path.isfile(table_fn):
                results_paths[analysis_type] = table_fn

        if not results_paths:
            logger.warning(
                f"No result tables found for {experiment_id} under "
                f"{os.path.join(str(results_dir), experiment_id, 'tables')}. "
                "It will contribute metadata but no statistics."
            )

        infos.append(AnalysisInfo(
            experiment_id=experiment_id,
            analysis_workbook=workbook,
            counts_path=counts_path,
            results_paths=results_paths,
        ))

    return infos


#### Experiment metadata ####


def doi_to_link(doi: str) -> str:
    """Format a DOI as a Markdown link. Does not resolve it."""
    if pd.isna(doi) or not doi:
        return ""
    for fragment in ("https", "http", "://", "doi.org/"):
        doi = doi.replace(fragment, "")
    return f"[{doi}](https://doi.org/{doi})"


def _find_year(reference: str) -> Union[int, str]:
    """Pull the publication year out of an MLA-style reference."""
    candidates = set()
    for part in reference.split(" (")[1:]:
        part = part.split(")")[0]
        try:
            year = int(part)
        except ValueError:
            continue
        if 1900 < year < 2100:
            candidates.add(year)

    if len(candidates) == 1:
        return candidates.pop()
    logger.warning(
        f"Could not settle on one year for reference {reference!r} "
        f"(candidates: {sorted(candidates)}). No citation created."
    )
    return "????"


def _short_citation(reference: str) -> str:
    """"{first_author} et al ({year})" from a full reference."""
    reference = normalise_text(reference)
    return f"{reference.split(',')[0]} et al ({_find_year(reference)})"


def tabulate_experiments_metadata(experiment_details: List[pd.Series]) -> pd.DataFrame:
    """One row per experiment, from each workbook's Experiment details sheet."""
    experiment_details = [
        details.copy().drop(np.nan, errors="ignore") for details in experiment_details
    ]
    for details in experiment_details:
        df_rename_columns(details, EXPERIMENT_COLUMN_RENAMER, inplace=True, axis="index")

    table = pd.DataFrame(experiment_details).reset_index(drop=True)

    # Older templates have no Reference field at all.
    if "Reference" not in table.columns:
        table.loc[:, "Reference"] = ""

    # Unpublished screens are cited by their experiment ID.
    unpublished = table.Reference.isna() | table.Reference.apply(
        lambda s: "internal" in str(s).lower()
    )
    table.loc[unpublished, "Reference"] = table.loc[unpublished].index
    table.loc[:, "Citation"] = ""
    table.loc[unpublished, "Citation"] = table.loc[unpublished].index
    table.loc[~unpublished, "Citation"] = (
        table.loc[~unpublished, "Reference"].apply(_short_citation)
    )

    unparsed = table.loc[table.Citation.str.contains("?", regex=False), "Experiment ID"]
    if len(unparsed):
        logger.warning(
            "Could not build a citation for these experiments:\n\t"
            + "\n\t".join(map(str, unparsed.values))
        )

    _deduplicate_citations(table)

    # DOI is not guaranteed to be present in every template revision.
    if "DOI" in table.columns:
        table["DOI"] = table["DOI"].apply(doi_to_link)
    else:
        logger.info("No DOI field in Experiment details; leaving it blank.")
        table["DOI"] = ""

    return table.infer_objects()


def _deduplicate_citations(table: pd.DataFrame) -> None:
    """Two screens from one paper become "Author et al (2020a)" and "(2020b)"."""
    citations = table["Citation"]
    if not citations.duplicated().any():
        return
    for citation, indices in citations.groupby(citations).groups.items():
        if len(indices) < 2:
            continue
        for n, index in enumerate(indices):
            table.loc[index, "Citation"] = citation.replace(")", chr(97 + n) + ")")


#### Comparison metadata ####


def get_treatment_str(sample_details: pd.DataFrame, ctrl: str, treat: str) -> str:
    """Describe in prose what changed between a control and a test sample."""
    treat_row, ctrl_row = sample_details.loc[treat], sample_details.loc[ctrl]

    if is_nt(treat_row.Treatment):
        chemical = ""
    elif is_nt(ctrl_row.Treatment) or ctrl_row.Treatment == treat_row.Treatment:
        chemical = treat_row.Treatment
    else:
        chemical = f"{ctrl_row.Treatment}{FAT_ARROW}{treat_row.Treatment}"

    # A knockout control against a wildtype treatment is assumed never to happen.
    if is_nt(treat_row.KO):
        knockout = ""
    elif is_nt(ctrl_row.KO) or ctrl_row.KO == treat_row.KO:
        knockout = f"{treat_row.KO}-KO"
        # A TP53 background is common enough not to be worth remarking on.
        if knockout == "TP53-KO":
            knockout = ""
    else:
        knockout = f"{ctrl_row.KO}-KO{FAT_ARROW}{treat_row.KO}-KO"

    if chemical and not knockout:
        return chemical
    if knockout and not chemical:
        return knockout
    if not knockout and not chemical:
        return "No treatment"
    if treat_row.Treatment == ctrl_row.Treatment:
        return f"{knockout} (with {chemical})"
    return f"{chemical} (in {knockout} cells)"


#: Sample fields copied onto a comparison, for the test sample and, prefixed with
#: "Control", for the control sample.
SAMPLE_FIELDS = (
    "Treatment", "KO", "Dose", "Growth inhibition %", "Days grown", "Cell line", "Notes",
)


def tabulate_comparisons(analysis_wb: AnalysisWorkbook) -> pd.DataFrame:
    """One row per comparison for a single experiment."""
    experiment_id = analysis_wb.expd["experiment_id"]
    logger.debug(f"Tabulating comparisons of {experiment_id}")

    rows = []
    for _, comparison in analysis_wb.control_groups.iterrows():
        ctrl = comparison["Control sample"]
        treat = comparison["Test sample"]

        info: Dict[str, object] = {"ControlSample": ctrl, "TestSample": treat}

        contrast = comparison.get("Contrast", "")
        if not contrast or len(str(contrast)) == 0:
            contrast = get_treatment_str(analysis_wb.wb["Sample details"], ctrl, treat)
        info["Contrast"] = contrast

        treat_row = analysis_wb.samples.loc[treat]
        ctrl_row = analysis_wb.samples.loc[ctrl]
        for field in SAMPLE_FIELDS:
            info[field] = treat_row.get(field)
            info["Control" + field] = ctrl_row.get(field)

        info["Experiment ID"] = experiment_id
        info["Library"] = analysis_wb.experiment_details.get("Library")
        info["Comparison ID"] = f"{experiment_id}.{ctrl}-{treat}"
        info["Timepoint"] = str(comparison["Group"]).split("_")[0]

        for field in ("Dose", "ControlDose"):
            if not pd.isna(info[field]):
                info[field] = str(info[field]).replace("uM", "μM")
        for field in ("KO", "ControlKO"):
            if pd.isna(info[field]) or info[field] in ("", "WT"):
                info[field] = "Wildtype"

        rows.append(info)

    comparisons = pd.DataFrame(rows)

    # "endpoint" is a common mis-entry for "endpoints".
    comparisons.loc[comparisons["Timepoint"] == "endpoint", "Timepoint"] = "endpoints"

    comparisons.loc[comparisons.Treatment.isna(), "Treatment"] = "No treatment"
    for column in ("Cell line", "Library", "Source"):
        if column not in comparisons.columns:
            comparisons.loc[:, column] = "Unspecified"
        comparisons.loc[comparisons[column].isna(), column] = "Unspecified"

    for internal, label in TIMEPOINT_LABELS.items():
        comparisons.loc[comparisons["Timepoint"] == internal, "Timepoint"] = label

    return comparisons


def create_metadata_tables(analysis_infos: List[AnalysisInfo]) -> MetadataTables:
    """Assemble the experiment and comparison metadata for every experiment."""
    comparisons = pd.concat(
        [tabulate_comparisons(info.analysis_workbook) for info in analysis_infos]
    )
    experiments = tabulate_experiments_metadata(
        [info.analysis_workbook.experiment_details for info in analysis_infos]
    )

    comparisons.set_index("Comparison ID", drop=False, inplace=True)
    experiments.set_index("Experiment ID", drop=False, inplace=True)

    return MetadataTables(comparisons=comparisons, experiments=experiments)


#### Statistics ####

#: Columns kept from each analysis type's results table.
_SHARED_STAT_COLUMNS = [
    "gene_id", "score", "fdr", "fdr10", "pos_p", "neg_p",
    "comparison_id", "analysis_type_id", "experiment_id",
]


def tabulate_statistics(info: AnalysisInfo) -> pd.DataFrame:
    """Flatten one experiment's result tables into StatTable rows."""
    tables = []

    for analysis_type, fn in info.results_paths.items():
        logger.debug(f"Tabulating {analysis_type.name} statistics from {fn}")
        try:
            stats = load_stats_csv(fn)
        except (pd.errors.ParserError, AttributeError):
            # An analysis that produced nothing leaves an unparseable table
            # behind; Chronos does this when no comparison qualified.
            logger.warning(f"Could not parse {fn}; skipping it.")
            continue

        stats.index.name = "gene_id"
        for comparison in stats.columns.levels[0]:
            table = stats[comparison].reset_index()
            df_rename_columns(
                table,
                {"lfc": "score", "normZ": "score", "fdr_log10": "fdr10"},
                inplace=True,
            )
            table.loc[:, "comparison_id"] = f"{info.experiment_id}.{comparison}"
            table.loc[:, "analysis_type_id"] = analysis_type.id
            table.loc[:, "experiment_id"] = info.experiment_id

            table = _select_stat_columns(table, analysis_type)
            if table is not None:
                tables.append(table)

    if not tables:
        return pd.DataFrame(columns=_SHARED_STAT_COLUMNS)

    return pd.concat(tables)


def _select_stat_columns(table: pd.DataFrame,
                         analysis_type: AnalysisType) -> Optional[pd.DataFrame]:
    """Fit one analysis type's columns to the shared statistics schema.

    Chronos reports a single score with no significance estimate, and manual
    analyses report one p-value. Both are broadcast across the significance
    columns so that the schema stays uniform.
    """
    if analysis_type.name in ("mageck", "drugz"):
        return table.loc[:, _SHARED_STAT_COLUMNS]

    if analysis_type.name == "chronos":
        source = "chronos_score"
    elif analysis_type.name == "manual":
        source = "pval"
    else:
        logger.warning(f"No column mapping for analysis type {analysis_type.name}.")
        return None

    if source not in table.columns:
        logger.warning(
            f"Expected a {source!r} column for a {analysis_type.name} table; "
            f"found {list(table.columns)}. Skipping."
        )
        return None

    if analysis_type.name == "manual":
        for column in ("fdr", "fdr10", "pos_p", "neg_p"):
            table[column] = table[source]
    else:
        for column in ("score", "fdr", "fdr10", "pos_p", "neg_p"):
            table[column] = table[source]

    return table.loc[:, _SHARED_STAT_COLUMNS]


#### Writing ####


def create_engine_with_schema(destination="sqlite://", echo=False) -> Engine:
    """Open a database and create any missing tables."""
    engine = create_engine(destination, echo=echo, connect_args={"timeout": 4000})
    TableBase.metadata.create_all(engine)
    return engine


def insert_records(table: Type[TableBase], records: List[dict],
                   session: Session) -> None:
    """Insert rows, dropping null values so that column defaults apply.

    Fails if a primary key already exists.
    """
    session.add_all([
        table(**{k: v for k, v in record.items() if not pd.isna(v)})
        for record in records
    ])


def get_gene_symbols_in_db(session: Session) -> set:
    """Every distinct gene symbol already recorded."""
    return {row[0] for row in session.query(GeneTable.symbol).distinct().all()}


def add_genes_from_symbols(symbols: Collection[str], organism: str,
                           session: Session) -> None:
    """Add placeholder gene records for symbols not yet in the database."""
    new_symbols = set(symbols).difference(get_gene_symbols_in_db(session))
    logger.debug(f"Adding {len(new_symbols)} new gene symbol(s)")
    insert_records(
        GeneTable,
        [dict(id=s, symbol=s, organism=organism) for s in new_symbols],
        session,
    )


def add_statistics(analysis_infos: List[AnalysisInfo], session: Session) -> None:
    """Insert every experiment's statistics into StatTable."""
    for info in analysis_infos:
        table = tabulate_statistics(info)
        if table.empty:
            logger.warning(f"No statistics to add for {info.experiment_id}.")
            continue

        organism = info.analysis_workbook.experiment_details.get("Organism")
        if pd.isna(organism) or not organism:
            logger.warning(
                f'{info.experiment_id} has no "Organism" in Experiment details; '
                "recording its genes as Human."
            )
            organism = "Human"

        add_genes_from_symbols(table.gene_id, organism, session)
        logger.info(f"Adding {table.shape[0]} statistics rows from {info.experiment_id}")
        insert_records(StatTable, table.to_dict(orient="records"), session)


def write_db_files(outdir: Union[str, Path], analysis_infos: List[AnalysisInfo],
                   metadata: MetadataTables, session: Session) -> None:
    """Add statistics, commit, then write the metadata CSVs."""
    add_statistics(analysis_infos, session)
    logger.info("Committing changes")
    session.commit()
    logger.info("Writing metadata tables")
    metadata.to_files(outdir)


def remove_experiments(metadata_tables: MetadataTables,
                       exp_ids: typing.Collection[str],
                       session: Session) -> MetadataTables:
    """Drop experiments from the statistics table and the metadata.

    Does not commit; the caller does that.
    """
    session.execute(
        sqlalchemy.delete(StatTable).where(StatTable.experiment_id.in_(list(exp_ids)))
    )
    return MetadataTables(
        comparisons=metadata_tables.comparisons.loc[
            ~metadata_tables.comparisons["Experiment ID"].isin(exp_ids)
        ],
        experiments=metadata_tables.experiments.loc[
            ~metadata_tables.experiments["Experiment ID"].isin(exp_ids)
        ],
    )


def confirm_deletion(existing: List[str]) -> None:
    """Get consent before replacing an existing database.

    Raises ConfirmationRequired rather than asking when there is nobody to ask.
    A non-interactive run must not silently destroy a database, but neither
    should it die on an EOFError from a prompt that could never be answered,
    which is what happens under `docker run` without `-it`.
    """
    print("These database files will be deleted:", file=sys.stderr)
    for fn in existing:
        print(f"    {fn}", file=sys.stderr)

    if not sys.stdin.isatty():
        raise ConfirmationRequired(
            "Refusing to replace an existing database without confirmation, and "
            "stdin is not a terminal so you cannot be asked. Re-run with "
            "--force-overwrite to replace these files, or allocate a terminal "
            "(`docker run -it ...`) to be prompted. To add to the existing "
            "database instead of replacing it, drop --new-db."
        )

    try:
        input("Press enter to continue, Ctrl+C to cancel. ")
    except (EOFError, KeyboardInterrupt):
        print(file=sys.stderr)
        raise ConfirmationRequired("Cancelled; the database was left alone.") from None


def create_database(outdir: Union[str, Path], analysis_infos: List[AnalysisInfo],
                    ask_before_deleting=True) -> None:
    """Build a database from scratch, replacing any already in outdir."""
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    existing = [
        fn for fn in glob.glob(str(outdir) + "/*") if Path(fn).name in DB_FILES
    ]
    if existing:
        if ask_before_deleting:
            confirm_deletion(existing)
        else:
            logger.info(f"Replacing {len(existing)} existing database file(s).")
        for fn in existing:
            os.remove(fn)

    engine = create_engine_with_schema(get_db_url(outdir))
    metadata = create_metadata_tables(analysis_infos)

    with Session(engine) as session:
        write_db_files(outdir, analysis_infos, metadata, session)

    logger.info(f"Created a database in {outdir}")


def update_database(db_dir: Union[str, Path], analysis_infos: List[AnalysisInfo],
                    update_experiments=False) -> None:
    """Add experiments to an existing database.

    Experiments already present are skipped, or replaced when
    update_experiments is set.
    """
    db_dir = Path(db_dir)
    missing = [fn for fn in DB_FILES if not (db_dir / fn).is_file()]
    if missing:
        raise FileNotFoundError(
            f"{db_dir} is not a database directory: missing {', '.join(missing)}. "
            "Use --new-db to create one."
        )

    metadata = MetadataTables.from_files(db_dir)
    engine = create_engine(get_db_url(db_dir))

    incoming = {info.experiment_id for info in analysis_infos}
    existing = set(metadata.experiments["Experiment ID"].unique())
    overlapping = incoming.intersection(existing)
    overlapping_str = ", ".join(sorted(overlapping)) or "none"

    with Session(engine) as session:
        if not update_experiments:
            logger.info(f"Already in the database, skipping: {overlapping_str}")
            analysis_infos = [
                info for info in analysis_infos
                if info.experiment_id not in overlapping
            ]
            modified = metadata
        else:
            logger.info(f"Already in the database, replacing: {overlapping_str}")
            modified = remove_experiments(metadata, list(overlapping), session)

        if not analysis_infos:
            logger.info("Nothing left to add.")
            return

        metadata = modified.join(create_metadata_tables(analysis_infos))
        write_db_files(db_dir, analysis_infos, metadata, session)

    logger.info(f"Updated the database in {db_dir}")


def remove_experiments_from_db(db_dir: Union[str, Path],
                               experiments_to_remove: Sequence[str]) -> None:
    """Remove experiments from a database on disk, committing the change."""
    db_dir = Path(db_dir)
    metadata = MetadataTables.from_files(db_dir)
    engine = create_engine(get_db_url(db_dir))

    with Session(engine) as session:
        modified = remove_experiments(metadata, experiments_to_remove, session)
        session.commit()

    modified.to_files(db_dir)
    logger.info(
        f"Removed {len(experiments_to_remove)} experiment(s) from {db_dir}"
    )
