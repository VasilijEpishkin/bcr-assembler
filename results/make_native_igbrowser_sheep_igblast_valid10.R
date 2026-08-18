library(igblastr)

ENV <- "/data/user/epishkin/conda/envs/bcr_env"

copy_browser_html <- function(out_html) {
  force(out_html)
  options(browser = function(url) {
    cat("igbrowser URL:", url, "\n")
    src <- sub("^file://", "", url)
    ok <- file.copy(src, out_html, overwrite = TRUE)
    cat("WROTE:", out_html, ok, "\n")
  })
}

is_present <- function(x) {
  !is.na(x) & x != "" & toupper(as.character(x)) != "NA"
}

compute_cdr3_coords <- function(df) {
  if (!"cdr3_start" %in% names(df)) df$cdr3_start <- NA_integer_
  if (!"cdr3_end" %in% names(df)) df$cdr3_end <- NA_integer_
  for (i in seq_len(nrow(df))) {
    if (is_present(df$cdr3_start[i]) && is_present(df$cdr3_end[i])) next
    cdr3 <- as.character(df$cdr3[i])
    seq <- as.character(df$sequence[i])
    if (!is_present(cdr3) || !is_present(seq)) next
    pos <- regexpr(cdr3, seq, fixed=TRUE)[1]
    if (!is.na(pos) && pos > 0) {
      df$cdr3_start[i] <- pos
      df$cdr3_end[i] <- pos + nchar(cdr3) - 1
    }
  }
  df
}

normalize_abstar_for_igbrowser <- function(df) {
  # abstar has AIRR-like calls plus FWR/CDR fields, but may lack cdr3_start/end.
  # Derive cdr3_start/end from the nucleotide cdr3 substring in sequence.
  df <- compute_cdr3_coords(df)
  if (!"sequence_id" %in% names(df) && "seq_id" %in% names(df)) df$sequence_id <- df$seq_id
  for (col in c("v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end", "cdr3_start", "cdr3_end")) {
    if (col %in% names(df)) df[[col]] <- suppressWarnings(as.integer(df[[col]]))
  }
  df
}

valid_for_igbrowser <- function(df) {
  required <- c("sequence", "sequence_id", "v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end", "cdr3_start", "cdr3_end")
  keep <- rep(TRUE, nrow(df))
  for (col in required) {
    if (!col %in% names(df)) return(rep(FALSE, nrow(df)))
    keep <- keep & is_present(df[[col]])
  }
  keep
}

make_native_browser <- function(in_tsv, out_html, n=10, mode=c("igblast", "abstar")) {
  mode <- match.arg(mode)
  copy_browser_html(out_html)
  x <- read.delim(in_tsv, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  if (mode == "abstar") x <- normalize_abstar_for_igbrowser(x)
  keep <- valid_for_igbrowser(x)
  y <- x[keep, , drop=FALSE]
  y <- y[1:min(n, nrow(y)), , drop=FALSE]
  cat("INPUT:", in_tsv, "\n")
  cat("valid rows:", sum(keep), "browser rows:", nrow(y), "\n")
  if (nrow(y) == 0) stop("No rows valid for igbrowser after filtering/conversion")
  print(y[, intersect(c("sequence_id", "productive", "v_call", "d_call", "j_call", "cdr3_aa", "v_sequence_start", "j_sequence_start", "cdr3_start", "cdr3_end"), colnames(y)), drop=FALSE])
  igbrowser(y)
}

make_native_browser(
  "/data/user/epishkin/results/PRJNA900592/annotator_compare/output/igblast/SRR22279249_200_igblast.tsv",
  "/data/user/epishkin/results/PRJNA900592/annotator_compare/igbrowser_SRR22279249_igblast_valid10.html",
  n=10,
  mode="igblast"
)
