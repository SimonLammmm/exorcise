#### The project version ####
#
# The single source of truth is the VERSION file at the root of the
# installation. Everything that reports a version reads it from there: both R
# entry points through here, the Python package through its __init__, and the
# image tag through docker/build-and-push.sh.
#
# This file deliberately depends on nothing, so that exorcise-trace can source
# it without pulling in the Bioconductor stack that the rest of R/ needs.


read_exorcise_version <- function(home) {
  path <- file.path(home, "VERSION")
  if (!file.exists(path)) {
    warning("No VERSION file at ", path, "; reporting the version as unknown.",
            call. = FALSE)
    return("unknown")
  }

  version <- trimws(readLines(path, n = 1L, warn = FALSE))
  if (length(version) == 0 || !nzchar(version)) {
    warning(path, " is empty; reporting the version as unknown.", call. = FALSE)
    return("unknown")
  }
  version
}
