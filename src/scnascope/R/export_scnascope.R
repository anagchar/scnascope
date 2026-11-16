# Write the small file scnascope opens, beside an ASCAT.sc result.
#
# Run it once at the end of an ASCAT.sc pipeline, where the result object is saved:
#
#   source(system("scnascope rscript", intern = TRUE))     # or the path of this file
#   save(res, file = "results/result_object.Rda")
#   export_scnascope(res, "results/result_object.Rda")
#
# This writes results/result_object.scnascope.h5ad (about 1-2% of the .Rda's size). scnascope
# looks for it first whenever it is given results/result_object.Rda (or the folder), so the
# large object never has to be read again:
#
#   scnascope --mode dna --tool ascatsc --input results/ --outdir viewer/
#
# Without the scnascope command where R runs, a compact results/result_object.scnascope.json.gz
# is written instead; scnascope opens that one too.

.scnascope_here <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(f)) {
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    f <- if (length(a)) sub("^--file=", "", a[1]) else NA
  }
  if (is.na(f)) getwd() else dirname(normalizePath(f))
})

export_scnascope <- function(res, path, scnascope = Sys.which("scnascope"), keep_json = FALSE) {
  if (!is.list(res) || is.null(res$allProfiles)) stop("res is not an ASCAT.sc result (no allProfiles)")
  # results/result_object.Rda -> results/result_object.scnascope.*
  stem <- sub("\\.(rds|rda|rdata)$", "", path, ignore.case = TRUE)
  json <- paste0(stem, ".scnascope.json")
  if (!exists("export_for_web", mode = "function"))
    source(file.path(.scnascope_here, "ascatsc_to_web.R"))
  app <- prepare_explorer_data(res)
  export_for_web(app, output_path = json)
  # gzip it: about 5x smaller, and scnascope reads .json.gz directly
  raw <- readBin(json, "raw", file.size(json))
  con <- gzfile(paste0(json, ".gz"), "wb", compression = 6)
  writeBin(raw, con)
  close(con)
  unlink(json)
  out <- paste0(json, ".gz")
  if (nzchar(scnascope)) {
    h5 <- paste0(stem, ".scnascope.h5ad")
    status <- system2(scnascope, c("convert", shQuote(out), "-o", shQuote(h5)))
    if (status == 0 && file.exists(h5)) {
      if (!keep_json) unlink(out)
      out <- h5
    } else {
      warning("scnascope convert failed; keeping ", out)
    }
  }
  message(sprintf("scnascope file: %s (%.1f MB)", out, file.size(out) / 1e6))
  invisible(out)
}

# the name before the rename from cnvscope
export_cnvscope <- export_scnascope
