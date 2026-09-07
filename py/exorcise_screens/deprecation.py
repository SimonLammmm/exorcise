"""Support for the command names this package used to be invoked by.

The tools were renamed when crispr_tools and crispr_screen_viewer were merged
into Exorcise. Every old name still works. Each prints a warning to stderr
naming its replacement, and is otherwise identical.

The mapping is declared here, in one place, so that the CLI, the shell wrappers
and the Docker entrypoint cannot drift apart.
"""

from __future__ import annotations

import sys
from typing import Callable, Dict, List, NamedTuple, Optional, Sequence, Tuple


class Deprecated(NamedTuple):
    """An old command name and what it became."""

    replacement: str
    #: Leading arguments the old form required which the new one does not take,
    #: e.g. `crispr-screen-viewer database` becomes `exorcise-database`.
    drop_leading: Tuple[str, ...] = ()
    #: Subcommands whose arguments need rearranging, not just stripping. Maps a
    #: subcommand to a function taking its arguments and returning the current
    #: equivalent.
    translate: Optional[Dict[str, Callable[[List[str]], List[str]]]] = None
    note: Optional[str] = None


def _translate_remove(args: List[str]) -> List[str]:
    """`remove DB_DIR EXP_ID...` becomes `--remove --out-dir DB_DIR EXP_ID...`.

    The old form took the database directory as its first positional argument;
    exorcise-database takes it as --out-dir.
    """
    if not args or args[0] in ("-h", "--help"):
        return ["--remove", "--help"]

    db_dir, targets = args[0], args[1:]
    if not targets:
        print(
            "WARNING: `crispr-screen-viewer remove` needs a database directory "
            "and then at least one experiment to remove.",
            file=sys.stderr,
        )
    return ["--remove", "--out-dir", db_dir, *targets]


# Keyed by the name the user typed.
DEPRECATED_COMMANDS: Dict[str, Deprecated] = {
    "ntByCycle": Deprecated("exorcise-trace"),
    "ntByCycle.R": Deprecated("exorcise-trace"),
    "count_reads": Deprecated("exorcise-count"),
    "count_reads.py": Deprecated("exorcise-count"),
    "crispr_pipeline": Deprecated("exorcise-analyse"),
    "crispr_pipeline.py": Deprecated("exorcise-analyse"),
    "crispr-screen-viewer": Deprecated(
        "exorcise-database",
        drop_leading=("database",),
        translate={"remove": _translate_remove},
        note=(
            "The `database` and `remove` subcommands are provided, as "
            "`exorcise-database` and `exorcise-database --remove`. The Dash web "
            "viewer, and the `launch`, `genes` and `test` subcommands, are no "
            "longer part of this package."
        ),
    ),
}

# The commands this package and its R sibling currently provide.
CURRENT_COMMANDS = {
    "exorcise": "Reannotate sequences against a genome and exome",
    "exorcise-trace": "Plot nucleotide frequency by sequencing cycle",
    "exorcise-count": "Count reads and map them to a guide library",
    "exorcise-analyse": "Run MAGeCK, DrugZ or Chronos over a screen",
    "exorcise-database": "Build a screen results database",
}


def warn_deprecated_name(invoked_as: str, stream=None) -> None:
    """Tell the user which name to use instead. Never raises."""
    entry = DEPRECATED_COMMANDS.get(invoked_as)
    if entry is None:
        return
    stream = sys.stderr if stream is None else stream

    print(
        f"WARNING: `{invoked_as}` is deprecated and will be removed in a future "
        f"release. It still works. Use `{entry.replacement}` instead.",
        file=stream,
    )
    if entry.note:
        print(f"         {entry.note}", file=stream)
    stream.flush()


def resolve_invocation(invoked_as: str,
                       args: Sequence[str]) -> Tuple[str, list]:
    """Map an old command name and its arguments onto the current equivalent.

    Returns (current_command_name, adjusted_args). Warns as a side effect when
    the name was a deprecated one. Unrecognised names are returned untouched, so
    this is safe to call unconditionally.
    """
    entry = DEPRECATED_COMMANDS.get(invoked_as)
    if entry is None:
        return invoked_as, list(args)

    warn_deprecated_name(invoked_as)
    args = list(args)

    # Subcommands whose arguments need rearranging rather than just stripping.
    if entry.translate and args and args[0] in entry.translate:
        subcommand = args.pop(0)
        return entry.replacement, entry.translate[subcommand](args)

    # `crispr-screen-viewer database ...` carried a subcommand that
    # `exorcise-database ...` does not.
    for expected in entry.drop_leading:
        if args and args[0] == expected:
            args.pop(0)
        elif not args or args[0] in ("-h", "--help"):
            # Let the real parser print its help.
            pass
        else:
            print(
                f"WARNING: `{invoked_as}` previously required the `{expected}` "
                f"subcommand. Continuing as `{entry.replacement} "
                f"{' '.join(args)}`.",
                file=sys.stderr,
            )

    return entry.replacement, args
