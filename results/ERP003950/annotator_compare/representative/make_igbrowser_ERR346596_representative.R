library(igblastr)

in_tsv <- "/data/user/epishkin/results/ERP003950/annotator_compare/output/igblast/ERR346596_200_igblast.tsv"
out_html <- "/data/user/epishkin/results/ERP003950/annotator_compare/representative/igbrowser_ERR346596_representative.html"
selected_ids <- c("ERR346596.2", "ERR346596.4", "ERR346596.20", "ERR346596.22", "ERR346596.28", "ERR346596.244", "ERR346596.941", "ERR346596.525", "ERR346596.903", "ERR346596.7", "ERR346596.27", "ERR346596.38", "ERR346596.44", "ERR346596.37", "ERR346596.45", "ERR346596.54", "ERR346596.56", "ERR346596.97", "ERR346596.68", "ERR346596.92")

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
