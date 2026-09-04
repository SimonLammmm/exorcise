"""SQLAlchemy schema for the screen results database.

Was crispr_screen_viewer/database.py, unchanged apart from comments and the
removal of unused imports.

Four tables are defined. The builder only writes to `gene` and `stat`;
comparison and experiment metadata goes into two gzipped CSVs alongside the
database rather than into SQL. `experiment` and `comparison` are still declared
because `stat` has foreign keys onto them, so the schema would not create
without them.
"""

from __future__ import annotations

from datetime import datetime
from typing import Optional

import sqlalchemy as sqla
from sqlalchemy import ForeignKey, orm
from sqlalchemy.orm import mapped_column

__all__ = ["TableBase", "GeneTable", "ExperimentTable", "ComparisonTable", "StatTable"]

mcol = mapped_column

# Declarative type shorthands.
MInt = orm.Mapped[int]
MFloat = orm.Mapped[float]
MStr = orm.Mapped[str]
# Nullable equivalents.
MIntN = orm.Mapped[Optional[int]]
MFloatN = orm.Mapped[Optional[float]]
MStrN = orm.Mapped[Optional[str]]


class TableBase(orm.DeclarativeBase):
    pass


class GeneTable(TableBase):
    __tablename__ = "gene"

    id: MStr = mcol(primary_key=True)
    symbol: MStr = mcol(unique=True)
    official_id: MStrN
    symbol_with_ids: MStrN
    organism: MStr


class ExperimentTable(TableBase):
    __tablename__ = "experiment"

    stringid: MStr = mcol(primary_key=True)
    date: orm.Mapped[Optional[datetime]]
    library: MStr
    doi: MStrN
    representation: MStrN
    moi: MStrN
    description: MStrN
    notes: MStrN
    reference: MStrN
    source: MStrN
    citation: MStrN


class ComparisonTable(TableBase):
    __tablename__ = "comparison"

    stringid: MStr = mcol(primary_key=True)
    experiment: MStr = mcol(ForeignKey(ExperimentTable.stringid))
    contrast: MStr
    treatment_label: MStr
    timepoint: MStr
    cell: MStr
    control_sample: MStr
    test_sample: MStr
    ko: MStr
    control_treatment: MStrN
    control_ko: MStrN
    dose: MStrN
    gi: MStrN
    days_grown: MIntN
    library: MStr
    notes: MStrN


class StatTable(TableBase):
    __tablename__ = "stat"
    __table_args__ = (
        sqla.UniqueConstraint("comparison_id", "gene_id", "analysis_type_id"),
    )

    # Composite primary key. The unique constraint above names these columns, so
    # renaming one means editing __table_args__ too.
    comparison_id: MStr = mcol(ForeignKey(ComparisonTable.stringid), primary_key=True)
    experiment_id: MStr = mcol(ForeignKey(ExperimentTable.stringid))
    gene_id: MStr = mcol(ForeignKey(GeneTable.id), primary_key=True)
    analysis_type_id: MInt = mcol(primary_key=True)

    score: MFloatN
    fdr: MFloatN
    fdr10: MFloatN
    pos_p: MFloatN
    neg_p: MFloatN
