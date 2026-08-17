# Project paths, resolved so scripts work whether the working directory is the
# repo root or src/ (RStudio tends to give one, `Rscript src/foo.R` the other).
#
# Every script starts with the same two-candidate bootstrap:
#   source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
# after which src_source() and the *_path() helpers are location-independent.

.find_root <- function() {
  for (p in c(".", "..", "../..")) {
    if (dir.exists(file.path(p, "data")) && dir.exists(file.path(p, "src")))
      return(normalizePath(p, winslash = "/"))
  }
  stop("cannot locate project root: expected sibling data/ and src/ directories")
}

PROJ_ROOT <- .find_root()

data_path <- function(...) file.path(PROJ_ROOT, "data", ...)
src_path  <- function(...) file.path(PROJ_ROOT, "src", ...)
out_path  <- function(...) {
  d <- file.path(PROJ_ROOT, "output")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  file.path(d, ...)
}

# Source a sibling script by bare name, from anywhere.
src_source <- function(file) source(src_path(file))
