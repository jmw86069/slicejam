# TODO items for Slicejam

## 17aug2023

* slicejam_analysis.Rmd

   * Make each step "robust" to error, with `tryCatch()` blocks,
   and displaying error messages, but continuing the workflow.
   * Enable custom statistical thresholds: MGM (already),
   `fold_cutoff`, `adjp_cutoff`.
   * Consider whether MASK regions should be excluded prior to
   statistical testing. Supporting case is that these regions may
   affect estimate of variance, and FDR adjustment (to minor degree imo).
   * DONE. Convert Sample Correlation heatmaps to `heatmap_se()`
   * DONE. Change all downstream sections to use iterate through
   multiple MGM thresholds (if present). I think this may make it much
   easier to test multiple MGM thresholds side-by-side, enabling us
   to test them in a more rigorous and systematic way.
   * limma contrasts now only include the "most appropriate" normalization,
   which starts with normalized data with blocking factor if present, then
   normalized data without blocking factor. To avoid batch-adjustment,
   imported data should not include column "Run".
   * PARTIAL. Debug extremely slow `limma contrasts` when using large number
   of rows, and with multiple MGM thresholds. Ultimately it is caused by
   using the `block` argument, instead of using batch-adjusted data.
   (Using batch-adjusted data is nearly equivalent except the model
   uses fewer degrees of freedom.)
   
      * It was caused by supplying `block` without `correlation`, and with
      too many rows. Internally `limma::lmFit()` calls
      `limma::duplicateCorrelation()` which scales "quadratically" (poorly)
      with increasing number of rows.
      * New argument `max_correlation_rows=10000` limits 10000 rows,
      which in testing produced the same `correlation` value.
      * Given `block` and `correlation` the `lmFit()` step took seconds,
      instead of 8-10 minutes.
      * Further, the `limma::voom()` step also called `lmFit()`, however
      it did not include `block` in the initial calculation of `weights`,
      as recommended by authors.
      * However authors now recommend a two-step approach:
      
         1a. call `voom()` without `block` to produce `weights`.
         1b. call `duplicateCorrelation()` with `weights` to produce
         `correlation` (using up to 10,000 rows).
         2a. call `voom()` with `block` and `correlation` to improve `weights`.
         2b. call `duplicateCorrelation()` with improved `weights` to improve
         `correlation` (using up to 10,000 rows).
      * Then given `block` and `correlation` the call to `lmFit()` should
      proceed as usual, but without the extreme delay caused by
      `duplicateCorrelation()`.

   * ACTIVE. Convert all sections to iterate consistently, so each type
   of plot has the same organization of tabs:
   
      1. MGM thresholds, which define the set of Peaks analyzed.
      2. Signal names (`assay_name`) when multiple normalizations are used.
      3. (optional) Contrast names, where individual contrasts are shown.
      4. (optional) Data centering, particularly for MA-plots and heatmaps.
   
   * DONE. Convert batch-adjusted contrasts to use non-batch-adjusted data,
   with batch as a blocking factor in `limma` with argument `block`.
   
      * Data must be stored making this change clear. It must therefore change
      `limma_batch_adjust_quantile_counts` to use `quantile_counts`,
      and add `blocking_values="Run"`.
      * New `"signal_label"` with `"quantile_counts (blocking by Run)"`.
      * All downstream steps must also display `signal_label`, while
      using `signal_name`.
      * Note heatmaps using `signal_name` should "ideally" use
      `centerby_colnames="Run"` or `"limma_batch_adjust_"` prefix,
      in order to visualize data consistent with the contrasts.
   
   * PARTIAL. Reduce data volume for `sestats` output. Some changes point
   upstream to `jamses::se_contrast_stats()`.
   * Enable or skip certain analysis steps:
   
      * Option to conduct contrasts using only normalized, or batch-adjusted data,
      saving time and data volume.
      * Consider option to skip peak annotation, which would also skip
      the peaks by region pie charts.
      * Consider option to save the `fc_se` object as RData for faster re-use.

   * Consider adding stacked bar graphs as an alternative to pie charts
   for "peaks by region". Pros and cons.
   

## 15aug2023

* PARTIAL. create helper function that creates a bash shell script that calls
the `run_slicejam.R` with proper environment.

   * `RHOME` should be defined using the current R session
   * `R_LIBS` should be defined using the current R user session
   * The file calling `Rscript` should therefore use the proper R
   executable, and the R packages installed for that user, hopefully
   making this script "portable" for use by others on the same linux system.
   * Need MS Windows shell script template equivalent.

* DONE. enable custom output directory from `slicejam_analysis.Rmd`

## 15mar2023

* Goal: Include Epilogos summary functions in the package.
Real goal: Allow use of chromatin state data.

   * Could be generalized to chromatin state data.
   * Consider R package that provides Epilogos summary data for re-use.
   * See raw data for hg19: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/epilogos/observed/
   * See also: https://epigenomegateway.wustl.edu/browser/?genome=hg19&datahub=https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/epilogos/imputed/qcat.json
   * Unclear where to obtain comparable data for hg38, mm39, other genomes,
   some liftOver has been done but the qcat.gz files are not available
   (or not found)
   * Perhaps hg38 here: http://compbio.mit.edu/epimap/
   * Annotate peaks using Epilogos data (only hg19 currently)

* Consider log-linear and/or chi-squared stats

   * enrichment, change in proportionality of peaks to feature_type_winner,
   or to epilogos state/state_group

## 09feb2023

RMarkdown updates:

0. Make output filenames and folders **smaller**. Please.

   * Option to provide custom output folder name.
   * Output HTML file should be:
   `"long_long_output_folder_analysis/slicejam_analysis.html"`, and not
   `"long_long_output_folder_analysis/long_long_output_folder.html"`

1. Remove the pre-MGM-filtered analysis steps. Instead, apply MGM filtering
to peaks before analysis, then run that analysis and follow-up.
2. Overall presentation:

   * consider some type of tab interface to present the different subsets
   and stat thresholds, e.g. adjP0.05, P0.05, MGM6 subset, etc.
   * add ability to specify which normalization to use:
   `"limma_batch_adjust_quantile_counts"`, or
   `"quantile_counts"` if batch adjustment should not be used.

2. Improve Venn diagrams

   * use `expand_fraction=0.2` so labels are visible.
   * clean up contrast names using `contrast2comp()`
   * remove prefix `"signal:"` and `"threshold:"` from figure titles.
   * attempt to add proportional Venn diagrams, scale down the `font_cex`.

3. Add option to pre-define group colors, so groups like
UL3, dH1A, and dH1B have consistently colored figures.
4. Bonus points for including epigenome regions as an alternative to
the `"feature_type"` Pie charts and bar charts.
5. Improve Pie charts:

   * clean up contrast names using `contrast2comp()`
   * split long contrast names across rows
   * remove "All Peaks" - since we have already filtered input peaks by MGM
   * assign colors to `"feature_type"` consistently across figures
   * order `"feature_type"` consistently in color legend

6. Improve "Annotate Broad Peaks"

   * consider renaming "Sliced Region Summary"
   * use `sestats_to_df()` format to report hit counts

7. Improve heatmaps

   

## 29nov2022

Functions to create command-line scripts with which to run the Rmarkown.

Constraints:

* Script must define RHOME based upon the active version of R and its path,
so that `Rscript` will call the appropriate R executable.
* Script must also define `.libPaths()` consistent with original user
settings to ensure the Rscript will run with the exact R library paths,
and with the exact version and executable of `R` itself.


## 25apr2022

* Migrate RMarkdown sections to use newer `jamses` functions

   * `design2colors()` to define colors by experiment grouping
   * `sestats_to_df()` to print `data.frame` summary of stat hits
   * `heatmap_se()`
   * `contrast2comp()` to shorten length of the contrast labels

* Reduce output filename length

   * enable custom output filename/subdirectory

* Improve clarity/consistency in section headings

   * MA-plots
   
      * panel ratios should be closer to 2:3 height:width, not 1:5
      * make labels shorter
   
   * Define clear "sestats types"
   
      * all peaks / adjP
      * all peaks / P-value
      * mgm peaks / adjP
      * mgm peaks / P-value
   
   * Volcano plot sections should match "sestats types"
   * Venn diagram sections should match "sestats types"
   
      * consider adding proportional Venn with `overlap_type="overlap"`
   
   * Pie charts of genomic regions
   
      * sections should match "sestats types"
      * subtitle should include number of hits represented
      * change "MGM-Filtered Peaks" to "All MGM Peaks"
   
   * Create bar chart format equivalent of pie chart data


## 31mar2022

* slicejam Rmarkdown

   * When the MA-plot layout is odd, for example 2 columns, 6 rows,
   the output plot dimensions should be adjusted to approximate
   square panels for each plot.
   
      * Affects: all MA-plots, Peak Width by Peak Signal

   * Contrasts:
   
      * need some way to limit contrasts to specific factors.
   
   * Sample correlations should probably not use hierarchical clustering
   * Section `"Limma Stats for All Peaks"` numbers do not appear to
   match the subsequent Volcano plots filtered to same criteria
   * Section `"Volcano for MGM Subset Peaks"` appears to be using the less
   stringent criteria, with `adjP=1`, instead of `adjP=0.05` as intended.
   
      * Correction: It displays Pval<0.05, then adjP<0.05 in adjacent plots,
      confusing. Should probably lead with adjP<0.05 for all plots, then
      a section with Pval<0.05 for all plots.
   
   * Venn diagrams:
   
      * contrasts should probably insert a line break between
      contrasts, otherwise labels get too wide to be displayed effectively.
      * Venn diagrams need a little table with total counts in each set.
      (Must be added to `venndir` package.)
      * Ideally place set labels outside, and all count labels inside,
      making sure that labels are near the top of the figure.
   
   * Pie charts:
   
      * Section `"Statistical Hits by Genome Region"`:
      Figure width is too narrow? Yikes, labels are not legible at all.
      Problem occurred with 2-column, 6-row layout.
      * Section `"MGM Subset Hits by Genome Region"` has great plot sizes
      for 2-column, 6-row layout.
      * These sections probably also need a stacked bar, or bars-beside
      format.
      * Contrast names probably need line break between group names,
      otherwise they are too wide.
   
   * Heatmap of Replicate Values
   
      * Needs to be wider per column/group
      * column_split labels need to be rotated 90 degrees, or they overlap.
      * Ideally exclude "Run" column when there is only one batch, no adjustment.
      * Might want to capture heatmap row order in the output tables so
      entries can be found in that table.
      * Add to vignette how someone could interrogate the heatmap row order
      as a drill-down to find which entries are being displayed.
   
   * Consider InteractiveComplexHeatmap?
   
      * Unclear if practical to embed a heatmap with 44k rows.
      * Main goal would be to allow zooming in to see rownames,
      gene assignment, genomic region, chromosome coordinates.
   
   * ATAC mode steps appeared to fail

* `mergeSplitCountBed.sh`

   * COMPLETE: allow flexibly providing gz BED files at input, 
   and with variable columns, by using only the first 3 BED columns
   * allow to define output filename directly, otherwise use previous approach
   * the featureCounts `DO_FC=1` and multiCovBed `DO_AUC=1` results are nearly
   identical, only differing by small amounts that suggest they handle the
   0-base/1-base BED entries differently; or they may handle unique/multimap
   reads differently. The majority of entries are identical integer values.
   
      * multiCovBed takes about 90 minutes, compared to about 1 minute for
      featureCounts - and produce mostly identical results. Not worth ever
      running multiCovBed
   
   * Eventually test process using bedgraph coverages with Adam Burkholder
   tool `make_heatmap` to see if true "area under the curve" might differ
   and be more accurate than featureCounts. If using Genrich ATAC coverage
   output (extend read ends to 100bp width) this coverage might be much more
   appropriate for analysis.


## 31jan2022


### design ideas:

* negative control lanes, e.g. Input, Vehicle, etc.

   * expect no real signal here
   * no normalization of these lanels
   * useful for identifying high noise, plot signal vs negative control
   and look for correlation where negative control has similar high signal,
   or high correlation
   * only allow comparisons within batch for negative control lanes
   * use the same MGM threshold, above which peaks potentially masked/filtered
   out of comparisons if the negative control signal is also above MGM

* peak calling Venn diagram across groups

   * requires peaks have signal above MGM in all samples of the group
   * total number of peaks analyzed
   * peaks passing MGM filter
   * peaks passing negative control filter
   * the current filtering only requires one replicate above MGM

* Re-run `slicejam_evaluation.Rmd` to confirm `iSamples` is properly used
at each step during batch adjustment.

* Add peak annotation to the `slicejam_evaluation.Rmd`?

   * alternative: Evaluation Mode for `slicejam_analysis.Rmd` that
   skips the statistical comparisons



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
