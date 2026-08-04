library(igblastr)

in_tsv <- "/data/user/epishkin/results/PRJNA900592/annotator_compare/output/igblast/SRR22279249_200_igblast.tsv"
out_html <- "/data/user/epishkin/results/PRJNA900592/annotator_compare/representative/igbrowser_SRR22279249_representative.html"
selected_ids <- c("SRR22279249.42", "SRR22279249.43", "SRR22279249.55", "SRR22279249.48", "SRR22279249.65", "SRR22279249.25", "SRR22279249.142", "SRR22279249.190", "SRR22279249.325", "SRR22279249.331", "SRR22279249.15", "SRR22279249.16", "SRR22279249.17", "SRR22279249.20", "SRR22279249.24", "SRR22279249.12", "SRR22279249.50", "SRR22279249.57", "SRR22279249.49", "SRR22279249.41")

options(browser = function(url) {
  cat("igbrowser URL:", url, "\n")
  src <- sub("^file://", "", url)
  ok <- file.copy(src, out_html, overwrite = TRUE)
  cat("WROTE:", out_html, ok, "\n")
})

x <- read.delim(in_tsv, sep="	", stringsAsFactors=FALSE, check.names=FALSE)
y <- x[x$sequence_id %in% selected_ids, ]
y <- y[match(selected_ids, y$sequence_id), ]
y <- y[!is.na(y$sequence_id), ]
required <- c("sequence", "sequence_id", "v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end", "cdr3_start", "cdr3_end")
keep <- rep(TRUE, nrow(y))
for (col in required) keep <- keep & !is.na(y[[col]]) & y[[col]] != ""
y <- y[keep, ]
cat("browser rows:", nrow(y), "\n")
print(y[, c("sequence_id", "v_call", "j_call", "cdr3_aa", "productive")])
igbrowser(y)
