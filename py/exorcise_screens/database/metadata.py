"""Analysis types, and the metadata tables that sit beside the database.

Extracted from crispr_screen_viewer/dataset.py. The DataSet class and the rest of
that module served the Dash application and are not included.
"""

from __future__ import annotations

import dataclasses
import typing
from pathlib import Path
from types import MappingProxyType
from typing import Collection, Union

import pandas as pd
from loguru import logger

COMPARISONS_CSV = "comparisons_metadata.csv.gz"
EXPERIMENTS_CSV = "experiments_metadata.csv.gz"
DB_FILENAME = "database.db"

#: Everything that makes up a database directory.
DB_FILES = (COMPARISONS_CSV, EXPERIMENTS_CSV, DB_FILENAME)


def get_db_url(db_dir) -> str:
    return f"sqlite:///{Path(db_dir)}/{DB_FILENAME}"


@dataclasses.dataclass
class MetadataTables:
    """The experiment and comparison metadata, held as two DataFrames."""

    comparisons: pd.DataFrame
    experiments: pd.DataFrame

    @classmethod
    def from_files(cls, directory: Union[Path, str]) -> "MetadataTables":
        directory = Path(directory)
        comparisons = pd.read_csv(directory / COMPARISONS_CSV)
        comparisons.set_index("Comparison ID", inplace=True, drop=False)
        experiments = pd.read_csv(directory / EXPERIMENTS_CSV)
        experiments.set_index("Experiment ID", inplace=True, drop=False)
        return cls(comparisons=comparisons, experiments=experiments)

    def to_files(self, directory: Union[Path, str]) -> None:
        directory = Path(directory)
        self.comparisons.to_csv(directory / COMPARISONS_CSV, index=False)
        self.experiments.to_csv(directory / EXPERIMENTS_CSV, index=False)

    def join(self, other: "MetadataTables") -> "MetadataTables":
        joined = MetadataTables(
            comparisons=pd.concat([self.comparisons, other.comparisons], axis="index"),
            experiments=pd.concat([self.experiments, other.experiments], axis="index"),
        )
        for attribute, id_column in (
            ("comparisons", "Comparison ID"),
            ("experiments", "Experiment ID"),
        ):
            table = getattr(joined, attribute)
            if table[id_column].duplicated().any():
                duplicated = table.loc[table[id_column].duplicated(), id_column]
                logger.warning(
                    f"Duplicate {id_column}(s) in {attribute} metadata after "
                    f"joining: {', '.join(map(str, duplicated.unique()))}"
                )
        return joined


@dataclasses.dataclass(frozen=True)
class AnalysisType:
    id: int
    name: str
    shortname: str
    label: str
    score_label: str


class _AnalysesTypes:
    """The available analysis types, addressable by name, short name or id.

    Use the ANALYSESTYPES singleton below rather than instantiating this.
    """

    def __init__(self, analyses_types: Collection[AnalysisType], default_type="drugz"):
        by_str = {v.name: v for v in analyses_types}
        by_str.update({v.shortname: v for v in analyses_types})
        self.by_str = MappingProxyType(by_str)
        self.by_id = MappingProxyType({v.id: v for v in analyses_types})
        self.list = tuple(analyses_types)
        self.default = self[default_type]

    def __getitem__(self, item: Union[str, int]) -> AnalysisType:
        if isinstance(item, str):
            return self.by_str[item]
        if isinstance(item, int):
            return self.by_id[item]
        raise ValueError(
            f"Analyses are addressed by str name or int id, not {type(item)}"
        )

    def __iter__(self) -> typing.Iterator[AnalysisType]:
        return iter(self.list)

    def __len__(self) -> int:
        return len(self.list)

    def str_to_id(self, name: str) -> int:
        return self.by_str[name].id


ANALYSESTYPES = _AnalysesTypes([
    AnalysisType(id=1, name="mageck", shortname="mag", label="MAGeCK",
                 score_label="Log2(FC)"),
    AnalysisType(id=2, name="drugz", shortname="drz", label="DrugZ",
                 score_label="NormZ"),
    AnalysisType(id=3, name="chronos", shortname="chr", label="Chronos",
                 score_label="Chronos score"),
    AnalysisType(id=4, name="manual", shortname="man", label="Manual",
                 score_label="Score"),
])
