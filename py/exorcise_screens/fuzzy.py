"""Mismatch-tolerant matching of read sequences onto a guide library.

This replaces `fuzzyTwoLists`, a Go program that was invoked through
`subprocess` and shipped as a committed binary. That binary was a macOS x86_64
build, so `--allow-mismatch` could not work inside the Linux container at all,
and it wrote three temporary files into the working directory without cleaning
them up.

The approach here avoids comparing every query against every library sequence.
Instead each library sequence is expanded into its own "neighbourhood": itself,
plus every sequence one substitution away. A query then matches by dictionary
lookup. That is O(L * k * 4) to build and O(Q) to query, rather than O(Q * L).

A query is only reported when it resolves to exactly one library sequence. If
one substitution could have produced it from two different guides, the read is
genuinely ambiguous and is discarded, which is the same rule the Go program
applied.
"""

from __future__ import annotations

from typing import Collection, Dict, Iterable, Set

from loguru import logger

BASES = ("A", "C", "G", "T", "N")

# Sentinel stored in the neighbourhood map for a sequence that more than one
# library guide could produce.
_AMBIGUOUS = object()


def _substitutions(sequence: str) -> Iterable[str]:
    """Every sequence exactly one substitution away from this one."""
    for position, original in enumerate(sequence):
        prefix = sequence[:position]
        suffix = sequence[position + 1:]
        for base in BASES:
            if base != original:
                yield prefix + base + suffix


def build_neighbourhood(lib_sequences: Collection[str],
                        max_mismatches: int = 1) -> Dict[str, object]:
    """Map every sequence within max_mismatches of a library sequence to it.

    Exact library sequences always win over near-matches, so a read that is
    itself a real guide is never treated as a typo of a different guide.
    """
    if max_mismatches not in (0, 1):
        raise ValueError(
            f"Only 0 or 1 mismatches are supported, got {max_mismatches}. "
            "Two or more mismatches would expand the library beyond what fits "
            "in memory and stops being specific enough to be useful."
        )

    exact: Set[str] = {str(s) for s in lib_sequences}
    neighbourhood: Dict[str, object] = {}

    if max_mismatches == 1:
        for sequence in exact:
            for variant in _substitutions(sequence):
                if variant in exact:
                    # Two library guides one substitution apart. The variant is
                    # a real guide, so it must not be claimed as a typo.
                    continue
                existing = neighbourhood.get(variant)
                if existing is None:
                    neighbourhood[variant] = sequence
                elif existing is not sequence and existing is not _AMBIGUOUS:
                    neighbourhood[variant] = _AMBIGUOUS

    # Exact matches are added last and unconditionally.
    for sequence in exact:
        neighbourhood[sequence] = sequence

    return neighbourhood


def map_allowing_mismatch(query_sequences: Collection[str],
                          lib_sequences: Collection[str],
                          max_mismatches: int = 1) -> Dict[str, str]:
    """Return {query: library_sequence} for queries that match unambiguously.

    Queries that match nothing, or that sit one substitution from two different
    library sequences, are omitted.
    """
    lib_sequences = [str(s) for s in lib_sequences]
    if not lib_sequences:
        logger.warning("The library contains no sequences; nothing can be matched.")
        return {}

    lengths = {len(s) for s in lib_sequences}
    if len(lengths) > 1:
        logger.info(
            f"Library guides vary in length ({min(lengths)}-{max(lengths)}); "
            "a read can only match a guide of its own length."
        )

    neighbourhood = build_neighbourhood(lib_sequences, max_mismatches)

    matched: Dict[str, str] = {}
    ambiguous = 0
    for query in query_sequences:
        query = str(query)
        target = neighbourhood.get(query)
        if target is None:
            continue
        if target is _AMBIGUOUS:
            ambiguous += 1
            continue
        matched[query] = target  # type: ignore[assignment]

    logger.info(
        f"Matched {len(matched)} of {len(query_sequences)} unmatched sequences "
        f"allowing {max_mismatches} mismatch(es)."
    )
    if ambiguous:
        logger.info(
            f"Discarded {ambiguous} sequence(s) that were one substitution from "
            "more than one library guide."
        )
    return matched
