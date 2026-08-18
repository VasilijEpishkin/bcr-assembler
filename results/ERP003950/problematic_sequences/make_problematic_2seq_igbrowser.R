# make_problematic_2seq_igbrowser.R
# Build igblastr::igbrowser HTML reports for the two problematic mouse sequences.
# Run inside OneQ/Jupyter task context, not through the SSH cgroup.

# OneQ tasks often have a non-writable /home/epishkin. Keep all R/igblastr
# caches on the shared writable /data/user volume before loading igblastr.
Sys.setenv(
  HOME = "/data/user/epishkin",
  XDG_CACHE_HOME = "/data/user/epishkin/.cache",
  R_USER_CACHE_DIR = "/data/user/epishkin/.cache/R"
)
options(igblastr_cache = "/data/user/epishkin/.cache/R/igblastr")
dir.create(getOption("igblastr_cache"), recursive = TRUE, showWarnings = FALSE)

.libPaths(c(
  "/data/user/epishkin/conda/envs/bcr_env/lib/R/library",
  "/data/user/epishkin/conda/envs/r45_igbrowser/lib/R/library",
  .libPaths()
))

suppressPackageStartupMessages(library(igblastr))

out_dir <- "/data/user/epishkin/results/ERP003950/problematic_sequences/igbrowser"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

igblast_tsv <- "/data/user/epishkin/results/ERP003950/problematic_sequences/output/igblast/problematic_mouse_2seq_igblast.tsv"
abstar_tsv  <- "/data/user/epishkin/results/ERP003950/problematic_sequences/output/abstar/problematic_mouse_2seq_abstar.tsv"

is_present <- function(x) {
  !is.na(x) & x != "" & toupper(as.character(x)) != "NA"
}

valid_for_igbrowser <- function(df) {
  required <- c(
    "sequence", "sequence_id",
    "v_sequence_start", "v_sequence_end",
    "j_sequence_start", "j_sequence_end",
    "cdr3_start", "cdr3_end"
  )
  keep <- rep(TRUE, nrow(df))
  for (col in required) {
    if (!col %in% names(df)) {
      cat("Missing required column:", col, "\n")
      return(rep(FALSE, nrow(df)))
    }
    keep <- keep & is_present(df[[col]])
  }
  keep
}

compute_cdr3_coords <- function(df) {
  if (!"cdr3_start" %in% names(df)) df$cdr3_start <- NA_integer_
  if (!"cdr3_end" %in% names(df)) df$cdr3_end <- NA_integer_
  for (i in seq_len(nrow(df))) {
    if (is_present(df$cdr3_start[i]) && is_present(df$cdr3_end[i])) next
    cdr3 <- as.character(df$cdr3[i])
    seq <- as.character(df$sequence[i])
    if (!is_present(cdr3) || !is_present(seq)) next
    pos <- regexpr(cdr3, seq, fixed = TRUE)[1]
    if (!is.na(pos) && pos > 0) {
      df$cdr3_start[i] <- pos
      df$cdr3_end[i] <- pos + nchar(cdr3) - 1
    }
  }
  df
}

normalize_abstar_for_igbrowser <- function(df) {
  # abstar TSV is AIRR-like but not native igblastr output.
  # Render only rows that contain enough coordinates for igbrowser.
  if (!"sequence_id" %in% names(df) && "seq_id" %in% names(df)) df$sequence_id <- df$seq_id
  if (!"v_call" %in% names(df) && "v_gene" %in% names(df)) df$v_call <- df$v_gene
  if (!"d_call" %in% names(df) && "d_gene" %in% names(df)) df$d_call <- df$d_gene
  if (!"j_call" %in% names(df) && "j_gene" %in% names(df)) df$j_call <- df$j_gene
  if (!"stop_codon" %in% names(df)) df$stop_codon <- NA
  if (!"vj_in_frame" %in% names(df)) df$vj_in_frame <- NA
  if (!"rev_comp" %in% names(df)) df$rev_comp <- NA
  df <- compute_cdr3_coords(df)
  for (col in c("v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end", "cdr3_start", "cdr3_end")) {
    if (col %in% names(df)) df[[col]] <- suppressWarnings(as.integer(df[[col]]))
  }
  # igbrowser expects numeric identity columns; abstar can leave D/J identity empty.
  # Use NA_real_ rather than character blanks so igbrowser passes type assertions.
  for (col in c("v_identity", "d_identity", "j_identity", "v_score", "d_score", "j_score", "v_support", "d_support", "j_support")) {
    if (!col %in% names(df)) df[[col]] <- NA_real_
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
  df
}

make_browser <- function(df, out_html, label) {
  force(out_html)
  options(browser = function(url) {
    cat("igbrowser URL:", url, "\n")
    src <- sub("^file://", "", url)
    ok <- file.copy(src, out_html, overwrite = TRUE)
    cat("WROTE:", out_html, ok, "\n")
  })
  cat("\n===", label, "===\n")
  cat("browser rows:", nrow(df), "\n")
  print(df[, intersect(c(
    "sequence_id", "productive", "stop_codon", "v_frameshift",
    "v_call", "d_call", "j_call", "cdr3_aa",
    "v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end", "cdr3_start", "cdr3_end"
  ), colnames(df)), drop = FALSE])
  igbrowser(df)
  if (!file.exists(out_html)) stop("igbrowser did not write HTML: ", out_html)
  cat("HTML_SIZE:", file.info(out_html)$size, "\n")
}

write_fallback_html <- function(df, out_html, label, reason) {
  esc <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  wanted <- intersect(c(
    "sequence_id", "productive", "stop_codon", "v_frameshift",
    "v_call", "d_call", "j_call", "v_identity", "j_identity", "v_support", "j_support",
    "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4",
    "fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa",
    "v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end",
    "fwr1_start", "fwr1_end", "cdr1_start", "cdr1_end", "fwr2_start", "fwr2_end",
    "cdr2_start", "cdr2_end", "fwr3_start", "fwr3_end", "cdr3_start", "cdr3_end",
    "fwr4_start", "fwr4_end"
  ), colnames(df))
  rows <- paste(vapply(wanted, function(col) {
    sprintf("<tr><th>%s</th><td><code>%s</code></td></tr>", esc(col), esc(df[[col]][1]))
  }, character(1)), collapse = "\n")
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>", esc(label), "</title>",
    "<style>body{font-family:Arial,sans-serif;line-height:1.35}table{border-collapse:collapse;width:100%}",
    "th,td{border:1px solid #ddd;padding:6px;vertical-align:top}th{width:180px;background:#f5f5f5;text-align:left}",
    "code{white-space:pre-wrap;word-break:break-all}</style></head><body>",
    "<h1>", esc(label), "</h1>",
    "<p><b>Fallback report, not native igbrowser rendering.</b></p>",
    "<p>Reason: ", esc(reason), "</p>",
    "<table>", rows, "</table></body></html>"
  )
  writeLines(html, out_html)
  cat("WROTE_FALLBACK:", out_html, file.info(out_html)$size, "\n")
}

# 1) Native IgBLAST/AIRR browser. Render each sequence separately because
# ko_seq_2 has no FWR4; igbrowser has an internal assertion failure on absent FWR4.
ig <- read.delim(igblast_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ig <- ig[ig$sequence_id %in% c("ko_seq_1", "ko_seq_2"), , drop = FALSE]
ig <- ig[match(c("ko_seq_1", "ko_seq_2"), ig$sequence_id), , drop = FALSE]
ig <- ig[!is.na(ig$sequence_id), , drop = FALSE]
ig_keep <- valid_for_igbrowser(ig)
cat("IgBLAST valid rows:", sum(ig_keep), "/", nrow(ig), "\n")

for (sid in c("ko_seq_1", "ko_seq_2")) {
  one <- ig[ig$sequence_id == sid, , drop = FALSE]
  if (nrow(one) == 0) next

  # Always write a simple field/table report so both sequences have the same
  # fallback-style representation for side-by-side inspection.
  fallback_html <- file.path(out_dir, paste0("problematic_", sid, "_igblast_fallback.html"))
  write_fallback_html(one, fallback_html, paste("IgBLAST", sid, "field report"), "explicit field/table report requested for comparable inspection")

  # Also try native igbrowser where the row shape is compatible.
  out_html <- file.path(out_dir, paste0("problematic_", sid, "_igblast_igbrowser.html"))
  if (!valid_for_igbrowser(one)[1]) {
    write_fallback_html(one, out_html, paste("IgBLAST", sid), "row lacks coordinates required by igbrowser")
    next
  }
  if ("fwr4" %in% names(one) && !is_present(one$fwr4[1])) {
    write_fallback_html(one, out_html, paste("IgBLAST", sid), "IgBLAST row has no FWR4; native igbrowser fails with .make_double_arrow_ascii(FWR4) assertion")
    next
  }
  tryCatch(
    make_browser(one, out_html, paste("IgBLAST", sid)),
    error = function(e) {
      cat("IGBROWSER_FAILED:", sid, conditionMessage(e), "\n")
      write_fallback_html(one, out_html, paste("IgBLAST", sid), conditionMessage(e))
    }
  )
}

# 2) abstar compatibility rendering only for ko_seq_1.
# ko_seq_2 has no abstar output row, so intentionally skip it.
ab <- read.delim(abstar_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ab <- normalize_abstar_for_igbrowser(ab)
ab <- ab[ab$sequence_id == "ko_seq_1", , drop = FALSE]
ab_keep <- valid_for_igbrowser(ab)
cat("abstar ko_seq_1 valid rows:", sum(ab_keep), "/", nrow(ab), "\n")
if (sum(ab_keep) > 0) {
  ab_one <- ab[ab_keep, , drop = FALSE]
  ab_out <- file.path(out_dir, "problematic_koseq1_abstar_igbrowser.html")
  tryCatch(
    make_browser(
      ab_one,
      ab_out,
      "abstar compatibility ko_seq_1 only"
    ),
    error = function(e) {
      cat("ABSTAR_IGBROWSER_FAILED:", conditionMessage(e), "\n")
      write_fallback_html(ab_one, ab_out, "abstar compatibility ko_seq_1 only", conditionMessage(e))
    }
  )
} else {
  cat("No abstar rows valid for igbrowser; skipped abstar HTML.\n")
}

cat("\nOUTPUT_DIR:", out_dir, "\n")
print(list.files(out_dir, pattern = "problematic_.*igbrowser.*html$", full.names = TRUE))
cat("DONE\n")
