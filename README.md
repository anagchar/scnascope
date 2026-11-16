# scnascope

Interactive viewer for single-cell copy-number calls, across modalities.

**Web app:** https://anagchar.github.io/scnascope: open a scnascope `.h5ad` or an ASCAT.sc
`.json(.gz)` export, or the demo dataset. Files are read in your browser and never uploaded.
The `scnascope` command (below) reads the tools' outputs and writes those `.h5ad` files.


| Modality | Tools |
|---|---|
| scDNA | ASCAT.sc, HiScanner, CHISEL, scphasedcn (long-read, phased), ASCENT |
| scATAC | epiAneufinder |
| scRNA | Numbat, InferCNV, CopyKAT |
| Multiome | Numbat-multiome (binned RNA + ATAC counts) → later paired per-modality calls |

Each tool's output is read into one canonical model (an AnnData object: cells × genomic bins,
one layer per value track) that records what each value *means* — integer copy number,
allele-specific copy number, gain/neutral/loss state, or a continuous ratio. The web UI picks
colour scales and views from that, not from the tool name, so adding a tool means adding a reader.

## Status

Everything the ASCAT.sc Explorer does, on the multi-tool model:

- upload screen (ASCAT.sc `ascat_data.json[.gz]`, `.h5ad`), demo dataset, Open / Panel / theme
- heatmap (WebGL): total or allele-specific (ASCAT.sc palette, its own dendrogram), hg38
  ideogram, drag-to-zoom with zoom chip, wheel pan, Ward.D2 dendrogram pruned to the zoom,
  annotation track, PNG/PDF export
- side panel: minimum segment size, heatmap height, BAF intervals, QC thresholds with
  "passing QC only", cell-type filter, cell search and list, ↑/↓ j/k, Esc, /, t
- cell profile with QC pills, BAF 5–95% ribbons, segment means and read depth; Lorenz curve;
  the cell's QC against all cells
- Cell QC tab (scatter explorer, Spearman matrix, distributions, statistics) and Dataset tab
- colours from the single-cell CNV literature: the SIGNALS / scgenome palette (Laks et al. 2019;
  Funnell et al. 2022); evenly graded, colour-blind-safe
  maps for ratios, BAF and counts (Crameri et al. 2020); Okabe–Ito for annotations
- "Figure" on the heatmap and the cell profile: publication figures with a live preview (journal
  widths, print theme, title, font size, dendrogram / ideogram / annotation / legend, whole genome
  or the zoomed region), as 300/600 dpi PNG or PDF

Readers (checked on real outputs): ASCAT.sc result objects (`.rds`/`.Rda`, via R) and web exports, HiScanner (`final_calls/`), CHISEL (`calls/calls.tsv` + `clones/mapping.tsv`). scphasedcn (a caller folder, or a sample folder holding `caller/` and `extract/`; checked on runs of its own simulator, not yet on real SPLONGGET data). scATAC: epiAneufinder (`epiAneufinder_results/results_table.tsv`; checked on tables written by its own R code, not yet on a real run).

## Install

```bash
conda env create -f environment.yml     # Python, Node.js (builds the interface), R for .rds/.Rda
conda activate scnascope
```

A conda package can be built from `conda/meta.yaml` (`conda build conda/`).

## Run

One command: the data type, the tool, its output folder, and where to write the results.

```bash
scnascope --mode dna --tool chisel    --input chisel_out/        --outdir results/
scnascope --mode dna --tool hiscanner --input hiscanner_output/  --outdir results/
scnascope --mode dna --tool ascatsc   --input ascatsc_results/   --outdir results/
scnascope --mode dna --chisel --input chisel_out/ --outdir results/      # shorthand
```

To read one specific file, for example an ASCAT.sc folder with several result objects:

```bash
scnascope --mode dna --tool ascatsc --input ascatsc_results/ --input-file results.rds --outdir results/
scnascope --mode dna --tool ascatsc --input-file /path/to/results.rds --outdir results/   # a full path alone
```

| `--mode` | `--tool` |
|---|---|
| `dna` | `ascatsc` (folder with the `.rds`/`.Rda` result or `ascat_data.json`), `hiscanner`, `chisel`, `scphasedcn` (the caller folder, the sample folder, or the `call_scn()` result object `<sample>_call_scn.rds`) |
| `rna` | `numbat` (its `out_dir`, or the folder holding `numbat_out/`), `infercnv` (its `out_dir`, or the folder holding `output/`), `copykat` (the folder with `*_copykat_CNA_results.txt`, or that file) |
| `atac` | `epianeufinder` (the run's output folder, holding `epiAneufinder_results/`) |
| `multiome` | `numbat` (a Numbat-multiome run, on genomic bins) |

It writes to `--outdir`, then opens the viewer:

- `<name>.h5ad`: the dataset; open it again in seconds with `scnascope serve results/<name>.h5ad`
- `<name>_cells.tsv`: one row per cell with its annotations (clone, cluster, cell type) and QC metrics
- `<name>_run.json`: what was run, on what, when, and what was written

Numbat: read at the final iteration (`clone_post_N`, `segs_consensus_N`, `geno_N`). Each cell
gets `p_cnv`, the probability of the most likely CNV segment over every bin (Numbat's `geno`),
and the segment's state as allele copies where that probability is at least 0.9 (del 1:0,
amp 2:1, loh 2:0, bamp 2:2, bdel 0:0; relative states, not absolute copy numbers), so it opens on
the allele-specific view. Bin-level runs (Numbat-multiome, "genes" named `chr1:0-1042457`) use
those bins and add the smoothed expression (`gexp_roll_wide`) as `expression`; gene-level runs average it per 1 Mb window, the genes placed with `bulk_clones_final.tsv.gz`. Cells are ordered by
Ward clustering of `p_cnv`, the matrix Numbat builds its phylogeny from; clone and tumour /
normal compartment are cell columns. A folder with several runs (e.g. `1MB/`, `220KB/`) is
refused with their names: point `--input` at one.

```bash
scnascope --mode multiome --numbat --input Numbat-multiome/1MB --outdir results/
```

InferCNV: read from its plain-text HMM results (no R needed): each cell is painted with its
subcluster's regions (`*.pred_cnv_regions.dat`, the Bayesian-filtered `HMM_CNV_predictions.*`
when present) on 1 Mb windows, HMM i6 states as copy number (1 → 0, 2 → 1, 3 → 2, 4 → 3, 5 → 4,
6 → 6); cell types from `cell_annotations.txt` and the subclusters are cell columns. With R on
the PATH, the smoothed expression InferCNV plots (`infercnv.png`) is added as `expression`, read
from `run.final.infercnv_obj` without the infercnv package, with the reference cells flagged.

Every scRNA / multiome tool (Numbat, InferCNV, CopyKAT) has this `expression` view next to its
calls: the tool's smoothed expression ratio, per window, as log2.

CopyKAT (the R package or its Python port): the windows' log ratios (`*_copykat_CNA_results.txt`)
as `expression`, and the prediction (`*_copykat_prediction.txt`: aneuploid / diploid; the Python
port's cluster and confidence too). A folder with several runs is refused with their names:
choose one with `--input-file`.

Expression-based ratios are small (CopyKAT ±0.15), so the heatmap's log-ratio colours span
± the 98th percentile of the data rather than ±1.

scphasedcn: normal mode paints each cell's reportable calls (arm, chromosome, 5-50 Mb tiles,
events) on 1 Mb windows over the neutral state, as allele-specific copy number (loss 1:0, gain
2:1, CN-LOH 2:0, ...; a male's X is haploid); tumour mode paints each assigned cell with its
clone's segments. With the extraction folder beside it (`extract/`, or the `dir` in the run
notes), per-cell read counts, the BAF (haplotype B / (A + B), shown under the profile in the allele-specific view) and a per-bin logR are added (each
window's read share against the null cells' median, or the reference pool's in tumour mode,
times its normal copy number, relative to the cell's mean copy number), so the cell profile
shows the bins around the calls as for the other tools. Its calls are calibrated
on null cells, so the viewer sets no default "Bins with CNA" filter for it.
When a run kept only the result object (`saveRDS(call_scn(...))`, no `write_scphasedcn()`),
point `--input-file` at the `.rds`, or `--input` at the sample folder with the `.rds` beside
it: R (with data.table; scphasedcn itself is not needed) writes its tables first, and
`<sample>/extract/` beside it is found for the read counts.

```bash
scnascope --mode dna --scphasedcn --input-file scphasedcn/COLO829BL_call_scn.rds --outdir results/
```

Cell annotations (cell types from RNA, ArchR clusters, samples, ...) are joined with
`--annotations cells.csv` (or `.tsv`): a barcode column (`cell`, `barcode`, or the first
column) and any others. Barcodes match with or without `cell-` prefixes, `-1` suffixes and
ArchR `sample#` prefixes; a `cell_type` column (or epiAneufinder's `annot`) becomes the
heatmap's annotation track.

```bash
scnascope --mode atac --epianeufinder --input epi_out/ --annotations cells.csv --outdir results/
```

`<name>` is the input folder's name (or the file's, for a lone `--input-file`), or `--name`.
Options: `--no-view` (only write the outputs), `--port` (the next free port is used if it is
taken, e.g. by another viewer), `--no-open`. `scnascope tools` lists the tools; `scnascope --help` shows everything.

Other commands: `scnascope serve [file or folder]` (no argument: the upload screen),
`scnascope convert <input> -o x.h5ad`, `scnascope demo`, `scnascope rscript`.

### ASCAT.sc: write the viewer file in the pipeline

Reading a large `.Rda`/`.rds` means starting R and loading the whole object, which takes
minutes on a network drive. Write a small viewer file beside it once, at the end of the
ASCAT.sc pipeline (where `scnascope` is installed, e.g. on the cluster):

```r
source(system("scnascope rscript", intern = TRUE))   # path of R/export_scnascope.R
save(res, file = "results/result_object.Rda")
export_scnascope(res, "results/result_object.Rda")    # -> results/result_object.scnascope.h5ad
```

Whenever scnascope is given `results/result_object.Rda` (or the `results/` folder), it first
looks beside it for `result_object.scnascope.h5ad` (or `result_object.h5ad`, or
`result_object.scnascope.json.gz`, which is written when the `scnascope` command is not
available to R) and opens that instead, in seconds. A viewer file older than the object is
ignored, and `--force` reads the object itself.

## Develop

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev]"
cd web && npm install
npm run dev                         # Vite on :5173, proxies /api to :8765
scnascope serve demo.h5ad --no-open  # in another shell
npm run build                       # writes the UI into src/scnascope/static
```

`pip install .` builds the interface itself when it isn't built yet and Node.js is available
(`hatch_build.py`).

Tests: `pytest`, `cd web && npm test`. Lint: `ruff check .`, `cd web && npm run lint`.

### The GitHub Pages site

The same interface, built to run without the Python server (`vite build --mode pages`):
`.h5ad` files are read with [h5wasm](https://github.com/usnistgov/h5wasm) and ASCAT.sc exports
in TypeScript (`web/src/local/`, with the QC metrics ported from `qc.py`; tests check both
against the Python). Deploy it to the `gh-pages` branch, with `scnascope` on the PATH (it writes
the demo dataset into the site):

```bash
cd web
npm run deploy        # build:pages (web/dist-pages + demo.h5ad), then push it to gh-pages
```

Preview it first with `npm run build:pages && python -m http.server -d dist-pages 8000`.

## Layout

- `src/scnascope/model.py`: the canonical dataset: modalities, value semantics, validation, .h5ad I/O
- `src/scnascope/genome.py`: chromosome sizes, bin grids, segments → bins
- `src/scnascope/analysis.py`: Ward.D2 / Manhattan clustering (total and allele-specific)
- `src/scnascope/qc.py`: per-cell QC metrics as ASCAT.sc's getQCs() computes them
- `src/scnascope/readers/`: one module per tool, plus a registry with format detection
- `src/scnascope/server.py`: FastAPI: dataset metadata and binary matrix chunks
- `src/scnascope/cli.py`: the command line (`--mode/--tool/--input/--outdir`, `serve`, `convert`)
- `web/`: Vite + React + TypeScript UI (WebGL2 heatmap)
