#### Common functions ####

# Create directory tree and write file
dirwrite <- function(x, dir, ...) {
  if (!file.exists(dir)) {
    dir.create(dirname(dir), recursive = T, showWarnings = F)
  }
  fwrite(x, dir, ...)
  invisible()
}

# Split strings
commasplit <- function(x) unlist(strsplit(x, ","))

# Return true if file doesn't exist or has zero size
file.not.exist.or.zero <- function(x) !file.exists(Sys.glob(paste0(x, "*"))[1]) | file.exists(Sys.glob(paste0(x, "*"))[1]) & file.size(Sys.glob(paste0(x, "*"))[1]) == 0

# fread with hanging glob
fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)