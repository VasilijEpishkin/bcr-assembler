#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: render_igbrowser.R input.airr.tsv output.html")
in_tsv <- normalizePath(args[[1]], mustWork = TRUE)
out_html <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  HOME = Sys.getenv("HOME", unset = "/tmp"),
  XDG_CACHE_HOME = file.path(Sys.getenv("HOME", unset = "/tmp"), ".cache"),
  R_USER_CACHE_DIR = file.path(Sys.getenv("HOME", unset = "/tmp"), ".cache", "R")
)
options(igblastr_cache = file.path(Sys.getenv("R_USER_CACHE_DIR"), "igblastr"))
dir.create(getOption("igblastr_cache"), recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(igblastr))

coordinate_cols <- c(
  "v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end",
  "j_sequence_start", "j_sequence_end", "cdr3_start", "cdr3_end",
  "fwr1_start", "fwr1_end", "cdr1_start", "cdr1_end", "fwr2_start", "fwr2_end",
  "cdr2_start", "cdr2_end", "fwr3_start", "fwr3_end", "fwr4_start", "fwr4_end"
)
numeric_cols <- c(
  "v_identity", "d_identity", "j_identity",
  "v_score", "d_score", "j_score",
  "v_support", "d_support", "j_support"
)

x <- read.delim(in_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
                quote = "", comment.char = "", na.strings = c("NA", "NaN"))
for (col in coordinate_cols) {
  if (!col %in% names(x)) x[[col]] <- NA_integer_
  x[[col]] <- suppressWarnings(as.integer(x[[col]]))
}
for (col in numeric_cols) {
  if (!col %in% names(x)) x[[col]] <- NA_real_
  x[[col]] <- suppressWarnings(as.numeric(x[[col]]))
}
for (col in c("sequence_id", "sequence", "locus", "v_call", "d_call", "j_call")) {
  if (!col %in% names(x)) x[[col]] <- ""
  x[[col]] <- as.character(x[[col]])
}
required <- c("sequence", "sequence_id", "v_sequence_start", "v_sequence_end",
              "j_sequence_start", "j_sequence_end", "cdr3_start", "cdr3_end",
              "fwr4", "fwr4_start", "fwr4_end")
stopifnot(all(required %in% names(x)))
keep <- rep(TRUE, nrow(x))
for (col in required) keep <- keep & !is.na(x[[col]]) & as.character(x[[col]]) != ""
x <- x[keep, , drop = FALSE]
if (nrow(x) == 0) stop("No representative rows valid for igbrowser")

copied <- FALSE
options(browser = function(url) {
  src <- sub("^file://", "", url)
  copied <<- file.copy(src, out_html, overwrite = TRUE)
  cat("igbrowser URL:", url, "\n")
  cat("WROTE:", out_html, copied, "\n")
})
cat("igblastr version:", as.character(packageVersion("igblastr")), "\n")
cat("input rows:", nrow(x), "\n")
print(table(x$report_library, useNA = "ifany"))
igbrowser(x)
if (!copied || !file.exists(out_html) || file.info(out_html)$size == 0) {
  stop("igbrowser did not create a non-empty durable HTML")
}
cat("output bytes:", file.info(out_html)$size, "\n")
provenance <- data.frame(
  tool = "igblastr::igbrowser",
  igblastr_version = as.character(packageVersion("igblastr")),
  R_version = R.version.string,
  input_rows = nrow(x),
  output_html = out_html,
  stringsAsFactors = FALSE
)
write.table(provenance, file.path(dirname(out_html), "igbrowser_provenance.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
