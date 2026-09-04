#!/usr/bin/env Rscript

# Deprecated name for exorcise-trace.R. Honoured, but warns.
#
# This is an R script rather than a shell wrapper so that `Rscript
# bin/ntByCycle.R` keeps working for anyone who invokes it that way.

message(
  "WARNING: `ntByCycle` is deprecated and will be removed in a future release. ",
  "It still works. Use `exorcise-trace` instead."
)

local({
  here <- dirname(normalizePath(
    sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]),
    mustWork = FALSE
  ))
  source(file.path(here, "exorcise-trace.R"))
})
