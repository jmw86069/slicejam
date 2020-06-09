# TODO items for Slicejam

## Enhancements

1. Change grouping input to "three-column style"

         Pattern            Batch         Group
         NOV14_p2w5_VEH     NOV14         p2w5_Veh
         NOV14_p4w4_VEH     NOV14         p4w4_Veh
         NOV14_UL3_VEH      NOV14         UL3_Veh
         NS644_UL3VEH       NS644         UL3_Veh
         NS50644_UL3VEH     NS50644       UL3_Veh
         NS644_p2w5VEH.     NS644         p2w5_Veh

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

4. Alternative data normalization method

   * Quantile (default)
   * median-scaling using `jamma::jammanorm()`
   * other convenient `limma` methods?

5. User-defined "active TSS" BED file?

   * substitutes for ATAC-mode annotation

6. Make size ranges for "Promoter region" and "TTS" user-tunable.

7. Export per-peak values

   * data matrix of raw, normalized, batch-adjusted values
   * data matrix of group mean values, raw, norm, batch
   * could be convenient place to add flags like "MGM peak", or "active promoter peak"

8. Export "active promoter TSS" as BED file

9. Improve pie chart labeling

   * reduce overlaps
   * make sure labels fit inside each plot panel


## R functions to port

* `volcanoPlot()`
* `statHitsAvenn()` - `vennPlotSets()`, `signedVennPlotSets()`
* `allNormStatsTests()`
* migrate from my `heatmap.3()` to `ComplexHeatmap::Heatmap()`


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
