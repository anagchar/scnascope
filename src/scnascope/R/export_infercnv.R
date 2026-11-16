# Export the smoothed expression of an InferCNV object (run.final.infercnv_obj) for scnascope.
#
#   Rscript export_infercnv.R run.final.infercnv_obj out_dir/
#
# The object is read as plain R attributes, so the infercnv package is not needed. Writes:
#   expr.f32     expr.data as little-endian float32, genes fastest (genes x cells, column-major)
#   genes.tsv    gene, chr, start, stop (gene_order, in the matrix's row order)
#   cells.txt    the matrix's column names
#   groups.tsv   cell, group, kind (reference / observation)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript export_infercnv.R run.final.infercnv_obj out_dir")
x <- readRDS(args[1])
out <- args[2]
e <- attr(x, "expr.data")
go <- attr(x, "gene_order")
if (is.null(e) || is.null(go)) stop("not an InferCNV object: no expr.data / gene_order")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

con <- file(file.path(out, "expr.f32"), "wb")
writeBin(as.numeric(e), con, size = 4, endian = "little")
close(con)
genes <- rownames(e)
go <- go[match(genes, rownames(go)), , drop = FALSE]
write.table(data.frame(gene = genes, chr = go$chr, start = go$start, stop = go$stop),
            file.path(out, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(colnames(e), file.path(out, "cells.txt"))

cells <- colnames(e)
grp <- function(lst, kind) {
  if (is.null(lst) || !length(lst)) return(NULL)
  do.call(rbind, lapply(names(lst), function(g) data.frame(cell = cells[lst[[g]]], group = g, kind = kind)))
}
groups <- rbind(grp(attr(x, "reference_grouped_cell_indices"), "reference"),
                grp(attr(x, "observation_grouped_cell_indices"), "observation"))
if (is.null(groups)) groups <- data.frame(cell = character(), group = character(), kind = character())
write.table(groups, file.path(out, "groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d genes x %d cells\n", nrow(e), ncol(e)))
