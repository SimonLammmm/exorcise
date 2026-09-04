"""Count reads in FASTQ/FASTA files and map them onto a guide library.

Merging assumes filenames of the form {sample}_L00?_R1_001.fastq[.gz] and that
the sequence of interest sits at the same offset in every read.

Two files are produced per sample: a dereplicated sequence count, and then a
single table of counts per guide once a library is supplied.

Was crispr_tools/count_reads.py.
"""

from __future__ import annotations

import datetime
import gzip
import os
from collections import Counter
from pathlib import Path, PosixPath, WindowsPath
from typing import Dict, Iterator, List, Sequence, Tuple, Union

import pandas as pd
from loguru import logger

from exorcise_screens.fuzzy import map_allowing_mismatch
from exorcise_screens.util import (
    ensure_parent_dir,
    file_exists_and_not_empty,
    is_temp_file,
    maybe_its_gz,
)

# Guides shorter than this are not counted: below about 16 bases a sequence stops
# being specific enough to identify a guide.
MIN_GUIDE_LENGTH = 16

FILE_FORMATS = {
    "a": (".fna", ".fasta", ".fa", ".fna.gz", ".fasta.gz", ".fa.gz"),
    "q": (".fastq", ".fastq.gz", ".fq", ".fq.gz"),
}
ALL_SUFFIXES = tuple(suffix for group in FILE_FORMATS.values() for suffix in group)


#### Reading ####


def read_fastq(file_obj) -> Iterator[str]:
    """Yield the sequence line of each FASTQ record, i.e. every fourth line."""
    for position, line in enumerate(file_obj):
        if position % 4 == 1:
            yield line.strip()


def read_fasta(file_obj, description_character=">") -> Iterator[str]:
    """Yield each FASTA record's sequence, joining wrapped lines."""
    sequence: List[str] = []
    for line in file_obj:
        if line.startswith(description_character):
            if sequence:
                yield "".join(sequence)
                sequence = []
        else:
            sequence.append(line.strip())
    if sequence:
        yield "".join(sequence)


def infer_file_format(fn: str) -> str:
    """'q' for FASTQ, 'a' for FASTA, by filename."""
    for fmt, suffixes in FILE_FORMATS.items():
        if fn.endswith(suffixes):
            return fmt
    raise ValueError(
        f"Cannot tell whether {fn} is FASTA or FASTQ from its name. "
        "Pass --file-type a or --file-type q."
    )


def count_reads_from_file(fn, window: Tuple[int, int], file_format="infer") -> Counter:
    """Count occurrences of the sequence found in `window` of each read.

    `window` is a (start, stop) pair of zero-based, stop-exclusive offsets.
    """
    if file_format == "infer":
        file_format = infer_file_format(fn)
    if file_format not in FILE_FORMATS:
        raise ValueError(f"Unknown file format {file_format!r}; expected 'a' or 'q'.")

    chop = slice(*window)
    counts: Counter = Counter()

    opener = gzip.open if fn.endswith(".gz") else open
    with opener(fn, "rt") as f:
        reader = read_fastq(f) if file_format == "q" else read_fasta(f)
        for sequence in reader:
            counts[sequence[chop]] += 1

    logger.info(
        f"{fn} window {tuple(window)}: {len(counts)} unique sequences, "
        f"{sum(counts.values())} reads"
    )
    return counts


def get_file_list(files_or_dirs) -> List[Path]:
    """Expand a mixture of files and directories into a list of files.

    Directories are listed one level deep, not recursively. Temporary and hidden
    files are skipped.
    """
    if isinstance(files_or_dirs, (str, PosixPath, WindowsPath)):
        files_or_dirs = [files_or_dirs]

    found: List[Path] = []
    for entry in (Path(f) for f in files_or_dirs):
        if entry.is_dir():
            found.extend(child for child in sorted(entry.iterdir()) if child.is_file())
        else:
            found.append(entry)

    return [fn for fn in found if not is_temp_file(fn)]


#### Counting ####


def guide_length_range(library_fn, seq_column: str) -> List[int]:
    """The shortest and longest guide in a library file."""
    library = pd.read_csv(maybe_its_gz(library_fn), sep=None, engine="python")
    if seq_column not in library.columns:
        raise RuntimeError(
            f'Library has no "{seq_column}" column. Columns are: '
            f"{list(library.columns)}. Set the sequence column with --seqhdr."
        )
    lengths = [len(str(s)) for s in library[seq_column]]
    return [min(lengths), max(lengths)]


def resolve_guide_lengths(library_fn, seq_column: str,
                          window: Sequence[int]) -> List[int]:
    """Work out which guide lengths to look for within the sliced region.

    A library whose guides are all the same length as the slice needs a single
    pass. Anything else needs a window slid along the slice, which is why this
    returns a range rather than one number.
    """
    slice_length = window[1] - window[0]

    if library_fn is None:
        return [slice_length, slice_length]

    lengths = guide_length_range(library_fn, seq_column)

    if lengths[0] < MIN_GUIDE_LENGTH:
        logger.warning(
            f"Library contains guides {lengths[0]} bases long. Not counting "
            f"anything shorter than {MIN_GUIDE_LENGTH} bases."
        )
        lengths[0] = MIN_GUIDE_LENGTH

    if lengths[0] != lengths[1]:
        logger.warning(
            f"Library guides vary in length ({lengths[0]}-{lengths[1]}). "
            "Sliding a window along the sliced region."
        )

    if lengths[0] > slice_length:
        raise RuntimeError(
            f"The shortest guide in the library is {lengths[0]} bases but the "
            f"slice is only {slice_length}. No guide can be found. Widen --slice."
        )

    if lengths[1] > slice_length:
        logger.warning(
            f"Slice length {slice_length} is shorter than the longest guide "
            f"({lengths[1]}). Guides longer than the slice will not be counted."
        )
        lengths[1] = slice_length
    elif lengths[1] < slice_length:
        logger.warning(
            f"Slice length {slice_length} exceeds the longest guide "
            f"({lengths[1]}). Sliding a window along the sliced region."
        )

    return lengths


def windows_within_slice(window: Sequence[int],
                         guide_lengths: Sequence[int]) -> List[Tuple[int, int]]:
    """Every sub-window of the slice that could hold a guide."""
    return [
        (start, start + length)
        for length in range(guide_lengths[0], guide_lengths[1] + 1)
        for start in range(window[0], window[1] - length + 1)
    ]


def group_files_by_sample(file_list: List[str], merge_samples: bool,
                          fn_split: str) -> Dict[str, List[str]]:
    """Map each output sample name to the input files that contribute to it."""
    if not merge_samples:
        return {fn.split(fn_split)[0].split("/")[-1]: [fn] for fn in file_list}

    # Lanes of one sample share everything up to _L001_.
    samples = {
        fn.split("_L001_")[0].split("/")[-1] for fn in file_list if "_L001_" in fn
    }
    grouped = {s: [fn for fn in file_list if s + "_L00" in fn] for s in samples}

    unmatched = set(file_list) - {fn for group in grouped.values() for fn in group}
    if unmatched:
        logger.warning(
            f"{len(unmatched)} file(s) have no _L001_ in their name and were not "
            f"merged into any sample: {', '.join(sorted(unmatched))}"
        )
    return grouped


def count_batch(files_or_dirs, window: Sequence[int], guide_lengths: Sequence[int],
                fn_prefix="", fn_suffix=".rawcount", fn_split="_R1_",
                merge_samples=False, file_type="infer",
                overwrite=False) -> List[str]:
    """Write one dereplicated sequence count file per sample.

    Skips any sample whose output already exists and is not empty, unless
    overwrite is set.
    """
    ensure_parent_dir(fn_prefix)

    file_list = [str(fn) for fn in get_file_list(files_or_dirs)]

    if file_type == "infer":
        ignored = [fn for fn in file_list if not fn.endswith(ALL_SUFFIXES)]
        if ignored:
            logger.warning(
                "These files do not look like FASTA/FASTQ and will be ignored:\n"
                + "\n".join(ignored)
            )
            file_list = [fn for fn in file_list if fn not in ignored]

    if not file_list:
        raise RuntimeError("No FASTA/FASTQ files to count.")

    grouped = group_files_by_sample(file_list, merge_samples, fn_split)
    sub_windows = windows_within_slice(window, guide_lengths)
    logger.info(
        f"Counting {len(file_list)} file(s) as {len(grouped)} sample(s), "
        f"guide lengths {guide_lengths[0]}-{guide_lengths[1]}, "
        f"{len(sub_windows)} window(s) per file."
    )

    written: List[str] = []
    for sample, sample_files in sorted(grouped.items()):
        outfn = _count_filename(fn_prefix, sample, fn_suffix)
        written.append(outfn)

        if file_exists_and_not_empty(outfn) and not overwrite:
            logger.info(f"Counts already exist, skipping: {outfn}")
            continue

        counts: Counter = Counter()
        for fn in sample_files:
            for sub_window in sub_windows:
                counts += count_reads_from_file(fn, sub_window, file_type)

        with open(outfn, "w") as f:
            for sequence, n in counts.most_common():
                f.write(f"{sequence}\t{n}\n")
        logger.info(f"Wrote {len(counts)} unique sequences to {outfn}")

    return written


def _count_filename(fn_prefix: str, sample: str, fn_suffix: str) -> str:
    if not fn_prefix:
        return f"{sample}.{fn_suffix}.txt"
    separator = "" if fn_prefix.endswith("/") else "."
    return f"{fn_prefix}{separator}{sample}{fn_suffix}.txt"


#### Mapping ####


def counts_table_from_files(file_list: List[Path], splitter=".raw",
                            remove_prefix=True) -> pd.DataFrame:
    """Assemble per-sample count files into one table indexed by sequence."""
    columns = {}
    for fn in file_list:
        fn = Path(fn)
        sample = fn.name.split(splitter)[0]
        if remove_prefix and "." in sample:
            sample = sample.split(".")[1]
        columns[sample] = pd.read_csv(fn, index_col=0, header=None, sep="\t")[1]

    logger.info(f"Sample columns: {list(columns)}")
    return pd.DataFrame(columns).fillna(0).astype(int)


def load_library(lib: Union[str, Path, pd.DataFrame], seqhdr: str, guidehdr: str,
                 genehdr: str) -> pd.DataFrame:
    """Read a guide library and check it has the three columns we need."""
    if isinstance(lib, (str, PosixPath, WindowsPath)):
        lib = str(lib)
        separator = "," if lib.endswith((".csv", ".csv.gz")) else "\t"
        lib = pd.read_csv(maybe_its_gz(lib), sep=separator)

    missing = [h for h in (seqhdr, guidehdr, genehdr) if h not in lib.columns]
    if missing:
        raise RuntimeError(
            f"Library is missing the column(s) {', '.join(missing)}.\n"
            f"Library columns: {list(lib.columns)}\n"
            f"Expected --seqhdr={seqhdr}, --guidehdr={guidehdr}, --genehdr={genehdr}"
        )

    if lib[guidehdr].duplicated().any():
        n = int(lib[guidehdr].duplicated().sum())
        logger.warning(
            f"{n} duplicated guide name(s) in the library. The first of each is "
            "kept and the rest discarded."
        )

    return lib


def map_counts(files_or_dirs, lib, seqhdr="seq", guidehdr="guide", genehdr="gene",
               report=False, splitter=".raw", remove_prefix=True, out_fn=None,
               allow_mismatch=False) -> pd.DataFrame:
    """Map counted sequences onto library guides.

    Returns a table indexed by guide name with a gene column and one column of
    counts per sample.
    """
    lib = load_library(lib, seqhdr, guidehdr, genehdr)
    lib = lib.set_index(seqhdr, drop=False)

    raw_counts = counts_table_from_files(
        get_file_list(files_or_dirs), splitter, remove_prefix
    )
    counts = raw_counts.reindex(lib[seqhdr], fill_value=0)

    if allow_mismatch:
        unmatched = raw_counts.loc[~raw_counts.index.isin(lib[seqhdr])].index
        now_matched = map_allowing_mismatch(unmatched, lib[seqhdr])

        if now_matched:
            gained = raw_counts.loc[list(now_matched)]
            gained.index = gained.index.map(now_matched)
            gained = gained.groupby(gained.index).sum()

            fraction = gained.sum().sum() / max(counts.sum().sum(), 1)
            logger.info(f"Mismatch matching recovered {fraction:.2%} more reads.")
            counts = counts + gained.reindex(counts.index, fill_value=0)

    # Guide names index the output because they are unique even when the same
    # sequence appears against two genes.
    counts.index = lib[guidehdr]

    if report:
        _report_mapping(counts, raw_counts)

    counts.insert(0, "gene", lib.set_index(guidehdr)[genehdr])
    counts = counts.drop_duplicates()

    if out_fn:
        ensure_parent_dir(out_fn)
        counts.to_csv(out_fn, sep="\t")
        logger.info(f"Wrote the counts table to {out_fn}")

    return counts


def _report_mapping(counts: pd.DataFrame, raw_counts: pd.DataFrame) -> None:
    total_raw = raw_counts.sum().sum()
    if not total_raw:
        logger.warning("No reads were counted at all.")
        return

    logger.info(f"{counts.sum().sum() / total_raw:.3%} of reads map to the library.")

    absent = int((counts == 0).all(axis=1).sum())
    logger.info(
        f"{absent / counts.shape[0]:.3%} ({absent}) of library guides were not found."
    )

    per_sample = (counts.sum() / raw_counts.sum()).apply(lambda n: f"{n:.3%}")
    logger.info(f"Mapping rate per sample:\n{per_sample.to_string()}")


def add_log_file(fn_prefix: str) -> None:
    """Tee the log into a timestamped file next to the outputs."""
    stamp = datetime.datetime.now().strftime("%y%m%dh%Hm%Ms%S")
    logger.add(f"{fn_prefix}.logfile.{stamp}.txt", level="INFO")
