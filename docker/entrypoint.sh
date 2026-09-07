#!/usr/bin/env bash
# Dispatch to one of the tools in the image.
#
# The R tools are invoked by absolute path because the R environment is not on
# the PATH. The Python tools are installed as console scripts in the Python
# environment, which is on the PATH, so they are called by name.
#
# The command names this image used to take still work; each warns and then runs
# its replacement. The mapping is the same as the one in
# py/exorcise_screens/deprecation.py.
set -euo pipefail

MAMBA_ROOT="${MAMBA_ROOT:-/opt/miniforge3}"
EXORCISE_HOME="${EXORCISE_HOME:-/exorcise}"
R_ENV="${MAMBA_ROOT}/envs/exorcise/bin"

usage() {
  cat <<EOF
Usage: docker run --rm -v .:/data <image> <command> [arguments...]

Commands:
  exorcise             Reannotate sequences against a genome and exome
  exorcise-trace       Plot nucleotide frequency by sequencing cycle
  exorcise-count       Count reads and map them to a guide library
  exorcise-analyse     Run MAGeCK, DrugZ or Chronos over a screen
  exorcise-database    Build or update a screen results database
  bash                 Open a shell in the container

Deprecated, still honoured:
  ntByCycle            -> exorcise-trace
  count_reads          -> exorcise-count
  crispr_pipeline      -> exorcise-analyse
  crispr-screen-viewer database
                       -> exorcise-database
  crispr-screen-viewer remove DB_DIR EXP_ID...
                       -> exorcise-database --remove --out-dir DB_DIR EXP_ID...

Run '<command> --help' for the arguments a command accepts.
EOF
}

deprecated() {
  echo "WARNING: \`$1\` is deprecated and will be removed in a future release. It still works. Use \`$2\` instead." >&2
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

command="$1"
shift

case "${command}" in
  exorcise)
    exec "${R_ENV}/Rscript" "${EXORCISE_HOME}/bin/exorcise.R" "$@"
    ;;
  exorcise-trace)
    exec "${R_ENV}/Rscript" "${EXORCISE_HOME}/bin/exorcise-trace.R" "$@"
    ;;
  exorcise-count)
    exec exorcise-count "$@"
    ;;
  exorcise-analyse)
    exec exorcise-analyse "$@"
    ;;
  exorcise-database)
    exec exorcise-database "$@"
    ;;

  # Deprecated names.
  ntByCycle | ntByCycle.R)
    deprecated "${command}" exorcise-trace
    exec "${R_ENV}/Rscript" "${EXORCISE_HOME}/bin/exorcise-trace.R" "$@"
    ;;
  count_reads | count_reads.py)
    deprecated "${command}" exorcise-count
    exec exorcise-count "$@"
    ;;
  crispr_pipeline | crispr_pipeline.py)
    deprecated "${command}" exorcise-analyse
    exec exorcise-analyse "$@"
    ;;
  crispr-screen-viewer)
    deprecated "crispr-screen-viewer" exorcise-database
    # This form took a subcommand. `database` and `remove` are provided; the
    # rest went with the Dash web viewer.
    case "${1:-}" in
      database)
        shift
        exec exorcise-database "$@"
        ;;
      remove)
        # The old form was `remove DB_DIR EXP_ID...`, with the database as the
        # first positional; exorcise-database takes it as --out-dir.
        shift
        if [[ $# -eq 0 || "$1" == "-h" || "$1" == "--help" ]]; then
          exec exorcise-database --remove --help
        fi
        db_dir="$1"
        shift
        exec exorcise-database --remove --out-dir "${db_dir}" "$@"
        ;;
      -h | --help | "")
        exec exorcise-database --help
        ;;
      launch | genes | test)
        echo "ERROR: the \`$1\` subcommand is no longer part of this image. Only \`database\` and \`remove\` are provided; the Dash web viewer was removed in 3.0.0." >&2
        exit 1
        ;;
      *)
        exec exorcise-database "$@"
        ;;
    esac
    ;;

  bash | sh)
    exec /bin/bash "$@"
    ;;
  -h | --help | help)
    usage
    exit 0
    ;;
  *)
    echo "Unknown command: ${command}" >&2
    echo >&2
    usage >&2
    exit 1
    ;;
esac
