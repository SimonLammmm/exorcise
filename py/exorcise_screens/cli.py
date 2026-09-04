"""Command line entry points.

Three commands are provided:

    exorcise-count       count reads and map them to a guide library
    exorcise-analyse     run MAGeCK, DrugZ or Chronos over a screen
    exorcise-database    build a screen results database

The names they used to have still work and are listed in deprecation.py.

Heavy imports are deliberately deferred into each command's body, so that
`--help` is immediate and `exorcise-count` never has to load TensorFlow.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional, Sequence

from exorcise_screens import __version__
from exorcise_screens.deprecation import resolve_invocation


def _add_version(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s (exorcise {__version__})"
    )


#### exorcise-count ####


def count_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="exorcise-count",
        description="Count unique sequences in FASTQ/FASTA files and, given a "
                    "library, map them onto guides. Filenames are assumed to look "
                    "like {sample}_L00?_R1_001.fastq[.gz].",
    )
    parser.add_argument(
        "files", nargs="+",
        help="FASTQ/FASTA files, or directories of them.",
    )
    parser.add_argument(
        "-s", "--slice", metavar="M,N", required=True,
        help="Region of each read holding the guide, as two zero-based "
             "comma-separated offsets, end exclusive. Required.",
    )
    parser.add_argument(
        "-p", "--prefix", metavar="PREFIX", required=True,
        help="Prefix for output files. May include a directory. Required.",
    )
    parser.add_argument(
        "--suffix", metavar="SUFFIX", default=".rawcount",
        help="Suffix for the per-sample count files, before .txt. [default: %(default)s]",
    )
    parser.add_argument(
        "--fn-split", metavar="STR", default="_R1_",
        help="Substring at which input filenames are split to derive the sample "
             "name. Ignored with --merge-samples. [default: %(default)s]",
    )
    parser.add_argument(
        "--merge-samples", action="store_true",
        help="Sum counts across files whose names differ only by lane.",
    )
    parser.add_argument(
        "--library", metavar="PATH", default=None,
        help="Guide library. When given, a mapped counts table is also written to "
             "<PREFIX>.counts.tsv.",
    )
    parser.add_argument(
        "-g", "--seqhdr", default="seq", metavar="COLUMN",
        help="Library column holding guide sequences. [default: %(default)s]",
    )
    parser.add_argument(
        "-j", "--guidehdr", default="guide", metavar="COLUMN",
        help="Library column holding guide names. [default: %(default)s]",
    )
    parser.add_argument(
        "-n", "--genehdr", default="gene", metavar="COLUMN",
        help="Library column holding gene names. [default: %(default)s]",
    )
    parser.add_argument(
        "--allow-mismatch", action="store_true",
        help="Also count reads one substitution from a library guide, when that "
             "guide is the only possible source.",
    )
    parser.add_argument(
        "--file-type", metavar="TYPE", choices=("a", "q", "infer"), default="infer",
        help="'a' for FASTA, 'q' for FASTQ, or infer from the filename. "
             "[default: %(default)s]",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Recount samples whose output already exists.",
    )
    parser.add_argument(
        "--no-log-file", action="store_true",
        help="Do not write a log file alongside the outputs.",
    )
    parser.add_argument(
        "--quiet", action="store_true", help="Only report warnings and errors.",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Report debug messages.",
    )
    # Accepted and ignored: counting no longer blocks for confirmation, so there
    # is nothing to skip.
    parser.add_argument("--just-go", action="store_true", help=argparse.SUPPRESS)
    _add_version(parser)
    return parser


def run_count(argv: Optional[Sequence[str]] = None) -> int:
    from loguru import logger

    from exorcise_screens.counting import add_log_file, count_batch, map_counts, \
        resolve_guide_lengths
    from exorcise_screens.util import set_loguru_level

    args = count_parser().parse_args(argv)

    level = "WARNING" if args.quiet else ("DEBUG" if args.debug else "INFO")
    set_loguru_level(logger, level)

    try:
        window = [int(n) for n in args.slice.split(",")]
        if len(window) != 2:
            raise ValueError
    except ValueError:
        print(
            f"--slice must be two comma-separated integers, got {args.slice!r} "
            "(for example --slice 0,20).",
            file=sys.stderr,
        )
        return 2
    if window[0] >= window[1]:
        print(
            f"--slice start must be less than its end, got {args.slice!r}.",
            file=sys.stderr,
        )
        return 2

    if args.library and not os.path.isfile(args.library):
        print(f"Library file not found: {args.library}", file=sys.stderr)
        return 2

    missing = [f for f in args.files if not os.path.exists(f)]
    if missing:
        print(f"Not found: {', '.join(missing)}", file=sys.stderr)
        return 2

    if not args.no_log_file:
        add_log_file(args.prefix)

    guide_lengths = resolve_guide_lengths(args.library, args.seqhdr, window)

    written = count_batch(
        files_or_dirs=args.files,
        window=window,
        guide_lengths=guide_lengths,
        fn_prefix=args.prefix,
        fn_suffix=args.suffix,
        fn_split=args.fn_split,
        merge_samples=args.merge_samples,
        file_type=args.file_type,
        overwrite=args.overwrite,
    )

    if args.library:
        map_counts(
            written, args.library,
            seqhdr=args.seqhdr, guidehdr=args.guidehdr, genehdr=args.genehdr,
            report=True, remove_prefix=True, splitter=args.suffix,
            out_fn=args.prefix + ".counts.tsv",
            allow_mismatch=args.allow_mismatch,
        )

    return 0


#### exorcise-analyse ####


def analyse_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="exorcise-analyse",
        description="Run CRISPR screen analyses described by an XLSX workbook or "
                    "a JSON configuration.",
    )
    parser.add_argument(
        "config_file", metavar="FILE",
        help="An .xlsx analysis workbook, or a .json configuration containing at "
             "least sample_reps, analyses and control_groups.",
    )
    parser.add_argument(
        "--counts", metavar="FILE/DIR", required=True, dest="counts_file",
        help="Counts file, or the directory holding the counts files the "
             "configuration names. Required.",
    )
    parser.add_argument(
        "--output-dir", metavar="PATH", default="./",
        help="Where the results directory is created. The experiment id and "
             "analysis version complete the path. [default: %(default)s]",
    )
    parser.add_argument(
        "--file-prefix", default="result",
        help="Prefix for every generated file. [default: %(default)s]",
    )
    parser.add_argument(
        "--skip-method", metavar="a,b", default=None,
        help="Comma-separated analysis methods to skip.",
    )
    parser.add_argument(
        "--run-groups", metavar="a,b", default=None,
        help="Comma-separated control groups to include. All by default.",
    )
    parser.add_argument(
        "--run-analyses", metavar="a,b", default=None,
        help="Comma-separated analysis names to include. All by default.",
    )
    parser.add_argument(
        "--analysis-version", default=None,
        help="Results are filed under this, if set.",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Rerun analyses whose output already exists.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Validate the configuration without running anything.",
    )
    parser.add_argument(
        "--dont-log", action="store_true", dest="dont_log", default=None,
        help="Do not write a log file.",
    )
    _add_version(parser)
    return parser


def run_analyse(argv: Optional[Sequence[str]] = None) -> int:
    from loguru import logger

    from exorcise_screens.pipeline import (
        PipelineOptionsError,
        load_configuration_file,
        process_arguments,
        run_analyses,
    )

    args = vars(analyse_parser().parse_args(argv))
    args.pop("version", None)

    # A directory given to --counts qualifies the filenames in the configuration;
    # a file overrides them entirely.
    counts_dir = ""
    counts_path = args.get("counts_file")
    if counts_path:
        if os.path.isdir(counts_path):
            counts_dir = counts_path
            del args["counts_file"]
        elif not os.path.isfile(counts_path):
            print(
                f"--counts {counts_path} is neither a file nor a directory.",
                file=sys.stderr,
            )
            return 2

    try:
        configuration = load_configuration_file(args["config_file"], counts_dir)
        # Anything given explicitly on the command line wins. --overwrite and
        # --dry-run default to False rather than None, so they always win; that
        # asymmetry is deliberate and long-standing.
        for key, value in args.items():
            if value is not None:
                configuration[key] = value

        run_analyses(**process_arguments(configuration))
    except PipelineOptionsError as error:
        logger.error(str(error))
        return 1

    return 0


#### exorcise-database ####


def database_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="exorcise-database",
        description="Build or update a screen results database from analysis "
                    "workbooks and the tables exorcise-analyse produced.",
    )
    parser.add_argument(
        "details_xlsx", nargs="+", metavar="DETAILS_XLSX",
        help="One or more analysis workbooks.",
    )
    parser.add_argument(
        "--out-dir", "-o", metavar="DIRECTORY", required=True,
        help="Where the database and metadata are written. Required.",
    )
    parser.add_argument(
        "--results-dir", "-r", metavar="DIRECTORY", default="results",
        help="Directory holding the analysis results. [default: %(default)s]",
    )
    parser.add_argument(
        "--counts-dir", "-c", metavar="DIRECTORY", default="counts",
        help="Directory holding the counts files. [default: %(default)s]",
    )
    parser.add_argument(
        "--filename-prefix", "-p", default="result",
        help="Prefix the analysis files were written with. [default: %(default)s]",
    )
    parser.add_argument(
        "--new-db", "-n", action="store_true",
        help="Create a new database. You are asked before existing database "
             "files are replaced, unless -f is also given.",
    )
    parser.add_argument(
        "--update-existing", "-u", action="store_true",
        help="Replace experiments already in the database instead of skipping them.",
    )
    parser.add_argument(
        "--force-overwrite", "-f", action="store_true",
        help="With --new-db, replace existing database files without asking.",
    )
    parser.add_argument(
        "--verbosity", "-v", metavar="N", type=int, default=1,
        help="0 for warnings only, 1 for info, 2 for debug. [default: %(default)s]",
    )
    _add_version(parser)
    return parser


def run_database(argv: Optional[Sequence[str]] = None) -> int:
    from loguru import logger

    from exorcise_screens.database.build import (
        ConfirmationRequired,
        create_database,
        get_paths,
        update_database,
    )
    from exorcise_screens.util import set_loguru_level

    args = database_parser().parse_args(argv)

    levels = ["WARNING", "INFO", "DEBUG"]
    verbosity = min(max(args.verbosity, 0), len(levels) - 1)
    if verbosity != args.verbosity:
        print(
            f"--verbosity must be 0, 1 or 2; clamping {args.verbosity} to "
            f"{verbosity}.",
            file=sys.stderr,
        )
    set_loguru_level(logger, levels[verbosity])

    missing = [fn for fn in args.details_xlsx if not os.path.isfile(fn)]
    if missing:
        print(f"Workbook(s) not found: {', '.join(missing)}", file=sys.stderr)
        return 2

    # A missing counts file or an unreadable workbook is a mistake in the
    # command, not a bug, so report it plainly rather than as a traceback.
    try:
        analysis_infos = get_paths(
            details_xlsx=args.details_xlsx,
            results_dir=args.results_dir,
            count_dir=args.counts_dir,
            analysis_filename_prefix=args.filename_prefix,
        )
    except (FileNotFoundError, RuntimeError, KeyError) as error:
        logger.error(str(error).strip("'"))
        return 2

    if not analysis_infos:
        print("No usable workbooks were given.", file=sys.stderr)
        return 2

    try:
        if args.new_db:
            create_database(
                args.out_dir, analysis_infos,
                ask_before_deleting=not args.force_overwrite,
            )
        else:
            update_database(
                args.out_dir, analysis_infos,
                update_experiments=args.update_existing,
            )
    except ConfirmationRequired as error:
        logger.error(str(error))
        return 1

    return 0


#### Dispatch ####

#: Every command this package provides, and the function that runs it.
COMMANDS = {
    "exorcise-count": run_count,
    "exorcise-analyse": run_analyse,
    "exorcise-database": run_database,
}


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Dispatch on the name the program was invoked as.

    This is what makes the old command names keep working: a wrapper named
    `count_reads.py` reaches the same code as `exorcise-count`, after a warning.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    invoked_as = os.path.basename(sys.argv[0]) if sys.argv else "exorcise-count"

    command, argv = resolve_invocation(invoked_as, argv)

    runner = COMMANDS.get(command)
    if runner is None:
        print(
            f"Unknown command: {command}. Available commands: "
            f"{', '.join(sorted(COMMANDS))}.",
            file=sys.stderr,
        )
        return 2

    return runner(argv)


#### Entry points ####
# Each of these is named in pyproject.toml. They exist so that a deprecated
# name resolves without depending on argv[0], which a shell wrapper or a
# container entrypoint may well have rewritten.


def _entry(command: str, argv: Optional[Sequence[str]], invoked_as: str) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    resolved, argv = resolve_invocation(invoked_as, argv)
    runner = COMMANDS[resolved if resolved in COMMANDS else command]
    try:
        return runner(argv)
    except KeyboardInterrupt:
        # Ctrl+C is a deliberate act, not a crash worth a traceback.
        print("\nInterrupted.", file=sys.stderr)
        return 130


def count_entry(argv: Optional[Sequence[str]] = None) -> int:
    return _entry("exorcise-count", argv, "exorcise-count")


def analyse_entry(argv: Optional[Sequence[str]] = None) -> int:
    return _entry("exorcise-analyse", argv, "exorcise-analyse")


def database_entry(argv: Optional[Sequence[str]] = None) -> int:
    return _entry("exorcise-database", argv, "exorcise-database")


def deprecated_count_reads(argv: Optional[Sequence[str]] = None) -> int:
    return _entry("exorcise-count", argv, "count_reads")


def deprecated_crispr_pipeline(argv: Optional[Sequence[str]] = None) -> int:
    return _entry("exorcise-analyse", argv, "crispr_pipeline")


def deprecated_crispr_screen_viewer(argv: Optional[Sequence[str]] = None) -> int:
    return _entry("exorcise-database", argv, "crispr-screen-viewer")


if __name__ == "__main__":
    sys.exit(main())
