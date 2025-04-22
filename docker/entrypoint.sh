#!/bin/bash
set -e

# Check the first argument to determine which application to run
case "$1" in
    exorcise)
        echo "Running exorcise:exorcise"
        exec /root/miniforge3/envs/exorcise/bin/Rscript /exorcise/bin/exorcise.R "${@:2}"
        ;;
    ntByCycle)
        echo "Running exorcise:ntByCycle"
        exec /root/miniforge3/envs/exorcise/bin/Rscript /exorcise/bin/ntByCycle.R "${@:2}"
        ;;
    count_reads)
        echo "Running crispr_tools:count_reads"
        exec /root/miniforge3/envs/crispr_tools/bin/count_reads.py "${@:2}"
        ;;
    crispr_pipeline)
        echo "Running crispr_tools:crispr_pipeline"
        exec /root/miniforge3/envs/crispr_tools/bin/crispr_pipeline.py "${@:2}"
        ;;
    crispr-screen-viewer)
        echo "Running crispr_screen_viewer:crispr-screen-viewer"
        exec /root/miniforge3/envs/crispr_screen_viewer/bin/crispr-screen-viewer "${@:2}"
        ;;
    *)
        echo "Unknown application: $1"
        echo "Usage: $0 {exorcise|ntByCycle|count_reads|crispr_pipeline|crispr-screen-viewer} [args...]"
        exit 1
        ;;
esac
