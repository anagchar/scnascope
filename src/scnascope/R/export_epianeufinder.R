# Export epiAneufinder's GC-corrected fragment counts (counts_gc_corrected.rds) for scnascope.
#
#   Rscript export_epianeufinder.R counts_gc_corrected.rds out_dir/
#
# The object is a data.table of windows x cells (the windows of results_table.tsv, in its order);
# it is read column by column, so neither data.table nor epiAneufinder is needed. Writes:
#   counts.f32   the counts as little-endian float32, windows fastest (one cell after another)
#   cells.txt    the column names (cell-<barcode>)
#   dims.txt     windows, cells

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript export_epianeufinder.R counts_gc_corrected.rds out_dir")
x <- readRDS(args[1])
if (!is.list(x) || !length(x)) stop("not an epiAneufinder counts table")
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
n <- length(x[[1]])
con <- file(file.path(out, "counts.f32"), "wb")
for (j in seq_along(x)) writeBin(as.numeric(x[[j]]), con, size = 4, endian = "little")
close(con)
writeLines(names(x), file.path(out, "cells.txt"))
writeLines(c(as.character(n), as.character(length(x))), file.path(out, "dims.txt"))
cat(sprintf("%d windows x %d cells\n", n, length(x)))
