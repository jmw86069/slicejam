# TODO items for Slicejam

## 25jan2022

Remove functions that were migrated into the `jamses` R package:

* se_normalize(), matrix_normalize(), groups_to_sedesign()
* se_contrast_stats(), run_limma_replicate(), etc.

Add `jamses` to package dependencies.


## 21jan2022: add epilogos annotations to genomic_regions_by_gtf

* Note it will only work for hg19 genome, when epilogos `GRangesList` is given.
* The epilogos data is 49 Mb, could be added as data to this package.
* Might require adding functions to parse epilogos data into `GRangesList`,
in the event someone found equivalent data from another genome.


## Bar chart alternative to pie charts

Pie charts were originally requested, but sometimes a stacked
bar chart could provide the same information in a more
compact and visually effective format.


## Minor convenience fixes

These are ideas without clear path to implement yet.

* `run_slicejam.R` should add file path when not specified, to protect
against directory changes.
* Helper function `update_slicejam()` which updates the R package
from github (if needed) then copies the two files from the package
directory to the current directory.


## Refactoring R code

* COMPLETE: Migrate from `heatmap.3()` to `ComplexHeatmap::Heatmap()`.
* COMPLETE: Move genomic regions to its own stand-alone R function.
* COMPLETE: Move volcano plot functions into this package (for now.)
* COMPLETE: Move Venn diagram functions to use `venndir::venndir()`.

   * Note some weirdness in font usage on linux systems, which
   appears to be mis-calculating font width, causing the name
   and directional label to overlap.


## Enhancements

1. COMPLETE: Change grouping input to "four-column style",
the first column (not shown) would have the filename, matching
the column header in the featureCounts `.fc` file.

         Sample             Run           Group
         NOV14_p2w5_VEH     NOV14         p2w5_Veh
         NOV14_p4w4_VEH     NOV14         p4w4_Veh
         NOV14_UL3_VEH      NOV14         UL3_Veh
         NS644_UL3VEH       NS644         UL3_Veh
         NS50644_UL3VEH     NS50644       UL3_Veh
         NS644_p2w5VEH      NS644         p2w5_Veh

1. In "mergeSplitCountBed" when run with "DRYRUN=1" it should
test for dependencies such as featureCounts.

1. Allow "Problem_Regions.bed" as optional input.

   * Peak slices are flagged when they overlap Problem_Regions
   * These peaks are not used during data normalization, stats tests

2. Annotate Broad Peaks

   * number of slices (already in place)
   * number of slices in "Problem_Regions.bed"
   * number of slices above MGM threshold
   * number of differential slices

3. Optionally filter for minimum peak slice width

   * peaks smaller than this width are not normalized nor tested

4.COMPLETE: Alternative data normalization methods

   * Quantile (default)
   * median-scaling using `jamma::jammanorm()`
   * median-scaling within sample group with `jamma::jammanorm()`
   * none - in case data is normalized outside slicejam
   * NOT COMPLETE: other convenient `limma` methods?

5.COMPLETE: User-defined "active TSS" BED file?

   * substitutes for ATAC-mode annotation
   * define `DETECTED_GENES` or `DETECTED_TX`

6. COMPLETE: Make size ranges for "Promoter region" and "TTS" user-tunable.

   * `UPSTREAM_PROMOTER`, `DOWNSTREAM_PROMOTER`, `UPSTREAM_TTS`, `DOWNSTREAM_TTS`

7. Export per-peak values

   * data matrix of raw, normalized, batch-adjusted values
   * data matrix of group mean values, raw, norm, batch
   * could be convenient place to add flags like "MGM peak", or "active promoter peak"

8. Export "active promoter TSS" as BED file

9. Improve pie chart labeling

   * reduce overlaps
   * make sure labels fit inside each plot panel


## R functions to port

* COMPLETE: `volcanoPlot()` - see `volcano_plot()`
* COMPLETE: `statHitsAvenn()` - `vennPlotSets()`, `signedVennPlotSets()`

   * see `venndir::venndir()`

* COMPLETE: `allNormStatsTests()`

   * see `slicejam::se_contrast_stats()`
   * see `slicejam::matrix_normalize()`, `slicejam::se_normalize()`

* COMPLETE: migrate from `heatmap.3()` to `ComplexHeatmap::Heatmap()`


## Things to test and demonstrate are handled in the workflow

* Compare differential peak results when using different peak callers

   * MACS2 peaks for a set of ChIP-seq data, run full pipeline, find hits
   * HOMER peaks for same ChIP-seq data, run full pipeline, find hits
   * Compare hits by direct region overlap
   * Compare hits by genes annotated

* Allow "Problem_Regions" to be included in peak slices

   * Show these extremely high count peaks in MA-plots
   * Show effect of normalization (do they adversely affect normalization?)
   * Show effect on sample correlation -- highest effect might be seen here
   * Show whether they become hits, if so, show how they can be filtered out

* Use "Problem_Regions.bed" to filter peak slices before normalization.

   * Show effect on normalization, compare to including Problem_Regions (above).
   * Show effect on sample correlation, should be much more accurate

* Show Broad Peak summary stats

   * Filter for Broad Peaks with >2 MGM slices
   * Filter for Broad Peaks with >2 differential MGM slices
   * Look for concordant and discordant pairs


## My notes from presentation 03jun2020

* "Takes you from peaks to annotated statistical hits"
* Better pie chart labeling to avoid overlaps
* Include mergeSplitCountBed

* Test differential workflow using public-downloaded peaks

   * all promoters
   * super-enhancers

* Needs way to define "Active genes", i.e. from ATAC-mode from
a separate pipeline.

## Comments/questions from presentation 03jun2020

* Compare different peak callers and see if it gives similar results.
* Use MACS2 with consistent sizes/parameters so they're about the same.
* Any thought to detecting problem peaks/regions in data? (Not yet.)
* Is the MGM threshold automatically set or user-defined? (User-defined.)
* Is promoter region pre-defined or numerically-defined? (Numeric)
And is it tunable command-line? (Yes)
* Asking if there are other ways to substitute the ATAC-mode aspect?
(Yes custom GTF)

* Trevor: How reliable is it to use series of histone marks to determine
active regions/promoters?

* Method to use RNA-seq data to determine "active transcription"

   * ATAC_NDR defines open promoters where there is no transcript,
   but where there could be induced transcription.
   * Jackson said of 13k active promoter ATAC sites, maybe only 500
   have no associated transcription
   
* Trevor: Are the promoter/TTS region sizes user-tunable? (Yes.)
