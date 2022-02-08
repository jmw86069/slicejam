# slicejam 0.0.37.900

## updates

* `slicejam_analysis.Rmd` was updated:

   * Section headings were slightly reorganized
   * to call `annotate_gr_by_genome_region()`
   instead of containing the logic inside the Rmarkdown file. The logic
   had been moved and tested outside Rmarkdown, since this function
   is also used by `genomic_regions_from_gtf`.
   * `feature_to_pie()` and `peak_pie_by_region()` were also
   moved outside the Rmarkdown file.
   * it now checks that the version of `slicejam` meets minimum requirements,
   in the unlikely event that the wrong `.libPaths()` location is used.
   * dependency on CBioRUtils was removed altogether. Woot!
   * `promoter_to_tss()` was moved out of `slicejam_analysis.Rmd` into
   a separate function.
  

## new functions

* `gg_pie_by_feature_type()` is a replacement pie chart function that
uses ggplot2 for improved aesthetics.
* `feature_to_pie()`, `peak_pie_by_region()` - draw pie charts by feature_type.
* `promoter_to_tss()` shrinks a promoter region down to its TSS position.
* `trim_multigenedesc()` is an internal function used to trim strings that
contain multiple gene descriptions, down to a maximum set number.

# slicejam 0.0.36.900

## updates

* Dependency was added to `jamses` into which several `SummarizedExperiments`
functions will be moved.
* `slicejam_analysis.Rmd` was updated to match the improved formatting
from `slicejam_eval.Rmd`. Eventually there will be one Rmd file with
switch that enables or disables the stats.
* Removed a bunch of functions that were migrated into the `jamses` package.

## possible bug fix

* `slicejam_analysis.Rmd` was updated to confirm the order of samples
in the batch adjustment step. The batch and group vectors are prepared
using `iSamples` but the `SummarizedExperiment` object is not ordered
using `iSamples`. In principle they should be the same order due to
upstream steps, but they should be defined in order to confirm.

# slicejam 0.0.35.900

## bug fixes

* Fixed typo bug in `genomic_regions_from_gtf` which caused user-defined
`outfile` not to be assigned, causing the script to end without exporting
results.

# slicejam 0.0.34.900

## bug fixes

* `matrix_normalize()` was failing on an edge case where numeric matrix
contained `NA` values, and `any(x <= floor)` failed because some
values were `NA` instead of `TRUE` or `FALSE`. The fix is
`any(x <= floor, na.rm=TRUE)` and again highlights one of the most
annoying R defaults: `na.rm=FALSE`.
It turns out `any()` will return `TRUE` if any value is `TRUE`, but
requires all values to be `FALSE` to return `FALSE` otherwise the
presence of one `NA` value will cause it to return `NA`.
* Updated several places in the code that checked numeric thresholds
that could possibly contain NA values.

# slicejam 0.0.33.900

## enhancements

* `printDebugList()` new argument `type="data.frame"` will return
a `data.frame` more suitable to Rmarkdown output.
* `slicejam_eval.Rmd` is a subset of `slicejam_analysis.Rmd` intended
to review only the peak slice count data, with no statistical comparisons.


# slicejam 0.0.32.900

## enhancements

* New process in `slicejam_analysis.Rmd` called `NormGroup` which
optionally allows specifying groups to use when `NORM="mediangroup"`.

   * This scenario is an edge case where a subset of samples are expected
   to have very low/absent signal, and therefore cannot be normalized
   relative to samples that have substantial signal. They are normalized
   within their own grouping, independent of other normalization groups.
   * By default, `NormGroup` is assigned by the sample group column in
   the `CURATION_TXT` file, however any column with `"NormGroup"` will
   be used instead if present.



# slicejam 0.0.31.900

## bug fixes

* Functions with `values()` were updated to proper `GenomicRanges::values()`.


# slicejam 0.0.30.900

## change in nomenclature

* Henceforth the term `"intergenic"` will be deprecated,
in its place we will use the term `"extragenic"` for these reasons:

   1. "intergenic" implies a region between genes, which would not
   represent regions between a gene and the end of a chromosome.
   2. "extragenic" refers to regions outside a gene body, regardless
   whether this region borders a gene or some other feature.
   3. "intergenic" may be confused with "intragenic" by others writing
   the term, which would mean the opposite of the intended definition.


## new functions

* `flatten_genome_regions()` takes the `genome_regions` data and
produces a single flat version with just the `feature_type` winners.


## changes to existing functions

* `annotate_gr_by_genome_region()` was updated to force consistent
output colnames, forcing `"gene_name"`, `"feature_type"`,
and when available in the source data: `"gene_id"`. Previously,
the code assumed `feature_type_colname="feature_type"` which would
have caused errors when the input colname was not `"feature_type"`.
* `annotate_gr_by_genome_region()` improved speed while checking
for gene_name values delimited by `","`.


# slicejam 0.0.29.900

## new object `SEDesign`

Object `SEDesign` are intended to contain the combination of:

* samples
* design matrix (samples x groups)
* contrasts matrix (groups x contrasts)

The intent is to allow subsetting samples, groups, and
to maintain integrity of the design and contrasts matrix.
For example, when one group in a contrast is removed,
the contrast should also be removed.

If samples are removed, and it results in removing a group,
it should also adjust the design and contrasts matrix data
accordingly.

* `validate_sedesign()` will update an existing `SEDesign` object,
optionally subset samples, groups, and contrasts.

Other helper functions

* `samples(SEDesign)` will return the sample identifiers
* `samples(SEDesign) <- ` will rename sample identifiers to values provided

* `groups(SEDesign)` will return the design groups
* `groups(SEDesign) <- ` will rename the design groups, including the groups in contrasts

* `design(SEDesign)` will return the design matrix
* `design(SEDesign) <- ` will reassign the design matrix

* `contrasts(SEDesign)` will return the contrast matrix
* `contrasts(SEDesign) <- ` will reassign the contrast matrix

## Next steps:

* include the `SEDesign` object as acceptable input to `se_contrast_stats()`
* change `se_contrast_stats()` to use `SEDesign` as input, and to convert
a la carte input (idesign, icontrasts, isamples) into `SEDesign`
to allow `SEDesign` methods to validate design and contrasts.



# slicejam 0.0.28.900

## changes to existing functions

* `ebayes2df()` was updated to handle contrast and interaction
contrast hit cutoffs more independently and correctly. Specifically,
the contrast and interaction contrast cutoff values are applied
independently.
* `run_limma_replicate()` arguments now reflect the cutoffs that
are passed down to `ebayes2df()`
* `se_contrast_stats()` now populates `hit_array` dimname `"Cutoffs"`
with distinct values for contrast cutoffs and interaction contrast
cutoffs. Previously they were handled in order but re-used labels
from only the contrast cutoffs.


# slicejam 0.0.27.900

## changes to existing functions

* `matrix_normalize()` logic for verbose output when
deciding whether to print `noise_floor_value`
was fixed to handle NA values.
* `se_normalize()` was updated to pass `verbose=verbose - 1` down to
`matrix_normalize()`.


# slicejam 0.0.26.900

## changes to existing functions

* `se_normalize()` now assigns the normalized matrix directly
in a way that retains any specific object attributes from
normalized matrix returned by `matrix_normalize()`.
This approach is useful for storing the `"nf"` normalization
factor values after `"jammanorm"` normalization for example.
Previously these attributes were lost during the assignment
step which used rownames,colnames to assign the subset
of samples and genes used during normalization.
* `se_normalize()` arguments now include the full set of
`params` recognized by `matrix_normalize()` in order to make
them visible in the function help text.
* `se_normalize()` and `matrix_normalize()` contain new
examples showing the workflow, with `farrisdata` read world
data used as a test set.


# slicejam 0.0.25.900

## changes to existing functions

* `volcano_plot()` changes:

   * changed so that it does not create a new blank page for each plot
   * adjusted `tophist=TRUE` logic so breaks are exactly aligned with
   volcano plot coordinates, and breaks always start exactly at the
   `fold_cutoff` edges. Internally it nowuses `graphics:::plot.histogram()`
   instead of `barplot()`. The x-axis tick marks use the same function
   `minorLogTicksAxis()` as the volcano plot, but without labels. The
   hope is that the tick mark alignment will make clear that the
   two axes are identical.

# slicejam 0.0.24.900

## changes to existing functions

* `volcano_plot()` changes:

   * `xlim` and `ylim` calculation with `range()` now includes
    `na.rm=TRUE` to account for NA values in input data.
    * `ylim` is now correctly calculated using `sig_floor`
    and `min_sig_range`.


# slicejam 0.0.23.900

## changes to dependencies

* `colorjam (>0.0.19.900)` to ensure `blend_colors()` is present.

## Normalization update

Normalization is configurable, initially adding median
normalization based upon `jamma::jammanorm()` which normalizes
the signal equivalent from an MA-plot, using median signal
at or above a minimum expression threshold. This method
is appropriate when normalization should not alter the
distribution of signal, instead should just shift the
signal up or down relative to other samples. The minimum
expression should be defined by the point in the MA-plot
where signal is clearly horizontal, and where the center
of such signal represents a stable frame of reference
for normalization.

* `run_slicejam.R` new argument `NORM` passed by environment
variable:

    * `NORM="quantile"` uses the default quantile normalization
    * `NORM="median"` uses median normalization
    * `NORM="mediangroup"` uses median normalization within group
    * `NORM="none"` performs no additional normalization
* `slicejam_analysis.Rmd` was refactored:

    * to apply the appropriate MGM to the initial unfiltered
    limma-voom analysis, where previously all peaks were tested
    but a default MGM=4 was applied to statistical hits.
    * to handle alternative
    normalizations, with appropriate plot labels indicating the
    method used.
    * to use `slicejam::volcano_plot()` instead of the proprietary
    `volcanoPlot()` used internally.
    * to use `ComplexHeatmap::Heatmap()` to migrate from using the
    never-released internal version of `heatmap.3()`.
    * New plot with the distribution of max group mean (MGM) values
    for putative hits, after filtering by P-value and fold change.
    This plot may be informative when selecting a reasonable MGM threshold.

## new functions

* `get_slicejam_args()` is a wrapper function to manage all
the input parameters sent by environment variable or directly
by R function argument. It is used for consistency in both
`run_slicejam.R` and `slicejam_analysis.Rmd`.

## mergeSplitCountBed.sh

* `mergeSplitCountBed.sh` was added to `exec/bin/mergeSplitCountBed.sh`

This file requires `bash` and is used as the command-line utility
to perform the peak merge-split-count workflow. In future this
workflow may be ported entirely into R, but currently it has
some utility being accessible on the command-line for use by
bench scientists.

Note that this script requires some other command-line utilities:

* bedtools - from Dr. Aaron Quinlan
* bedops
* featureCounts - from the Subread tool suite
* samtools
* perl - no specific modules are required

# slicejam 0.0.22.900

## updates to Rmarkdown

* In `slicejam_analysis.Rmd` the section `add_gene_names`
was updated for more streamlined logic
when adding gene description to an existing row.
* The file `run_slicejam.R` was updated to print the `ATAC`
variable, and to copy `slicejam_analysis.Rmd` to the current
working directory if it does not already exist.

# slicejam 0.0.21.900

## updates to existing functions

* `matrix_normalize()` new argument `enforce_norm_floor=TRUE`
re-applies the `floor` value to the normalized output data.
The effect is that when input data is zero, the output data
is also zero and is not allowed to be adjusted above or below
zero. The assumption for this default is that a value of zero
is not a measurement but represents the lack of a measurement.
Similarly, the intent of `floor` is a numeric threshold at or
below there is no confidence in the reported measurement, therefore
values at or below this threshold are treated as equivalent
to the threshold for the purpose of downstream analyses.
* `matrix_normalize()` help docs were edited with substantially
more text, and set of examples with simulated data showing
several characteristics of normalization and MA-plots using
`jamma::jammaplot()`.


# slicejam 0.0.20.900

## bug fixes

* `volcano_plot()` fixed bug where applying fold change ceiling
forced points to exponentiated max values instead of log2
values. This bug affected visual display and not the
determination of points above the fold ceiling threshold,
which caused points to be placed outside the viewing area.

## updates to existing functions

* `volcano_plot()`

    * Help docs were updated to describe all function arguments.
    * New argument `mar_min` to enforce at least a minimum margin
    around the plot panel, to allow room for block arrows.
    * New arguments `min_sig_range` and `min_fold_range` so the
    plot will at least span a reasonable range of significance and
    fold changes, even with data that has low significance
    or small fold changes.


# slicejam 0.0.19.900

## new functions

Three new functions associated with peak overlap calculations:

* `peakoverlap_calcs()` performs all-by-all overlap calculations,
using several relevant metrics. The example includes preparing
a multiple-panel heatmap to compare metric.
* `peakoverlap_heatmap()` is a helper function that prepares
a heatmap using one metric from `peakoverlap_calcs()`. It
can optionally include a label in each cell with other metrics.
* `cell_fun_label()` is a helper function that adds a text 
label to each cell in a heatmap via `ComplexHeatmap::Heatmap()`.

## removed defunct files

* `"slicejam_setup.md"` and `"slicejam_setup.html"` were removed as
the information is defunct. The replacement is `"how_to_slicejam.md"`.
* Marked a number of R functions "defunct" with description comment
that points to the recommended R function replacement. In most
cases these functions were used internally so the changes are
not visible to R users.


# slicejam 0.0.18.900

## updates

* `genomic_regions_from_gtf` was updated to handle `genes` and/or
`tx` input as files.
* All Venn functions will be replaced with new R package `venndir`.
* `genomic_regions_from_gtf()` calls `GenomicRanges::reduce()` on the
TTS regions per gene, which effectively collapses all overlapping
TTS regions for a gene into one contiguous region. This change
reduces the number of TTS regions to 1/3 the original size, while not
reducing content.
* `genomic_regions_from_gtf()` now only keeps unique TSS positions
used to calculate promoter regions for each gene. It still keeps
each unique promoter region distinct (different from TTS regions)
because the slicejam pipeline expects to be able to determine
the TSS from the promoter region during the "ATAC mode" annotation,
where a TSS whose promoter region contains an ATAC peak is
considered an "active TSS".


## bug fixes

* `genomic_regions_from_gtf()` was updated to fix several bugs
when supplied with `detectedGenes` and `detectedTx`.
* `genomic_regions_from_gtf()` does not try to save `RData` when
supplied with `detectedGenes` or `detectedTx`. The `RData` is
intended to represent the full GTF without any subset operations.
* `genomic_regions_from_gtf()` fixed incomplete TTS range, it only
correctly extended upstream the TTS and not downstream the TTS.
* `genomic_regions_from_gtf()` now uses `GenomicRanges::flank()`
instead of `GenomicRanges::promoters()` to extend the promoter region,
to keep the subset genes/transcripts consistent when creating
promoter and TTS regions.
* `slicejam_analysis.Rmd` now properly uses custom `upstream_promoter`
and `downstream_promoter` during the "ATAC mode" annotation step,
previously it used defaults 2000, and 200, respectively.


# slicejam 0.0.17.900

## updates

* updated `genomic_regions_from_gtf()` after command-line testing. More
work forthcoming to make it as smooth a standalone tool as possible.

# slicejam 0.0.16.900

## bug fixes

* `slicejam_analysis.Rmd` had one remnant reference to colname
`"Group"` that needed to be changed to `group_colname`, in
a lesser-used MA-plot section (`ma_q_group_run`) that organizes
panels by groups and batches when there are multiple replicates per
batch. This layout is intended to show whether the batches
are balanced across groups, and the relative distributions
relative to each batch.

# slicejam 0.0.15.900

## changes to existing functions

* `genomic_regions_from_gtf()`

    * returns the RData file
    saved, as an attribute `"rdata_file"`
    * new argument `rdata_file`
    allows defining a specific file to save the resulting
    `genome_regions` GRanges object. When that file exists,
    the data will be loaded and used as-is.
    * new argument `force_refresh` will force re-processing
    the data even if the output `rdata_file` exists
    * new argument `mask_regions` will add a column `"mask"`
    with `TRUE` or `FALSE` indicating whether each `GRange`
    overlaps any mask region.
    * new argument `save_bed` will save the genome_regions
    in BED format.

# slicejam 0.0.14.900

## new functions

* `annotate_gr_by_genome_region()` annotates GRanges using
`genome_regions` as returned by `genomic_regions_from_gtf()`.
This function makes it more feasible to use a command-line
function to annotate a BED file with a GTF file.

## venn functions in progress

Several functions related to Venn diagrams are in development,
partial forms are included but not ready for full release yet.

# slicejam 0.0.13.900

## changes to existing functions

* `volcanoPlot()` was renamed to `volcano_plot()`
* `volcano_plot()` received many updates to argument names,
removing most of the camelCase, and reducing overall arguments.
(More reduction to be done.)
* `blockArrowMargin()` now adjusts itself appropriately when the
figure aspect ratio is not 1:1 -- which is dramatic when the ratio
is 1:2 or 2:1.

## new function

* `find_colname()` is a utility function, probably moving to jamba,
which finds the best suitable matching colname in a `data.frame` given
one or more text patterns. Specifically here it is good for finding
P-value columns, or adjusted P-value columns; fold change or log2fold
columns, etc.

# slicejam 0.0.12.900

## new functions

The Venn update! (Not fully functional yet, functions being migrated.)

* `signed_venn_counts()` assembles various forms of directional
Venn counts, colors, and UTF8 arrows as labels.
* `list2im_opt()`, `list2im_signed()` are optimized forms of
`multienrichjam::list2im()` -- which as it turns out is painfully
slow even using `arules` `transactions` operations, when the input
list contains more than 10,000 entries each. Operations were taking
20+ seconds, now complete in less than 1 second. Could probably
still be faster.
* `match.matrix()` fill an unmet need, borrowing code from `unique.matrix()`
* `match.im()` is an optimized match for incidence matrix data,
converting them into integer scores per row -- but obviously limiting
the input incidence matrices to about 16 columns. Sufficient for this
purpose anyway.
* `blend_colors()` - finally ported a custom function, and after much
testing feel good about how it operates. It emulates the effect of
blending paint colors -- subtractive color mixing, and using a
red-yellow-blue color wheel instead of red-green-blue. Most driving
use cases work as intended now:

   * blue + yellow = green
   * red + yellow = orange
   * red + blue = purple
   * red + yellow + blue = grayish brown (as expected from paints)
   
   Note that it also handles color weights, which is used when some
   colors have partial transparency. It also handles blending more than
   two colors -- as shown in the examples above.
* `mean_angle()` takes one or more angles in degrees, and returns
the average angle, based upon the average unit vector. It also
handles weights, so certain angle unit vectors can be scaled
accordingly. (This function is used by `blend_colors()` to produce
an average color hue.) It also returns the radius, which works for
color hue blending, because the radius represents how much color
chroma (saturation) remains after taking the average.
* `symbol2utf8()` is a wrapper function to return one of various
UTF-8 arrow characters: `"upArrow"`, `"downArrow"` and several more.

# slicejam 0.0.11.900

## new functions

* `volcanoPlot()` and supporting function `blockArrowMargin()`,
`gradient_rect()`, `logAxis()`. Will be migrated probably to
the `"jamma"` package, for now is here to migrate away from
custom functions not stored in an R package.

# slicejam 0.0.10.900

## new functions

* `ebayes2dfs()` converts the limma eBayes fit into a list of
`data.frame` which also include a column for hits, based upon
one or more statistical thresholds.
* `log2fold_to_fold()` and `fold_to_log2fold()` convert to and
from signed fold change, where fold change -4 represents log2 fold
change -2.
* `mark_stat_hits()` takes a `data.frame` of stats results, and marks
statistical hits by three types of criteria: detected, changing, and
significant.

# slicejam 0.0.9.900

## enhancements

* `slicejam_analysis.Rmd` allows more custom settings:

    * upstream_promoter,downstream_promoter,upstream_tts,downstream_tts
    * `detectedTx` defined by environment variable `DETECTED_TX`,
    used to subset the GTF transcripts used when defining genome_regions.
    * `detectedGenes` defined by environment variable `DETECTED_GENES`,
    used to subset the GTF transcripts used when defining genome_regions.

## new functions

* `se_contrast_stats()` is a workhorse function, wrapper to run
multiple statistical contrasts on a `SummarizedExperiment` object,
primarily using `limma` methods.
* `run_limma_replicate()` is a wrapper to run `limma` contrasts on
a single data matrix, applying one or more contrasts in bulk.
* `voom_jam()` is a Jam-specific minor modification to `limma::voom()`
for the specific scenario where there may be NA values, which causes
the `stats::lowess()` function to return erroneous results.
* `handle_na_values()` handles a numeric matrix that contains `NA`
values, in a group context. Three options: `"full"` replaces only
`NA` values when the entire group is `NA`; `"partial"` replaces `NA`
only when the group contains some non-`NA` values. This function is
aimed at platform technologies where missing values are `NA` and
have particular meaning to that platform. Notably, when all values
of a group are `NA`, limma will return `NA` for all related contrasts.
When any values are not `NA` it is sometimes preferable to keep the
value without adding any replacement. When an entire group is `NA`
is can be useful to replace `NA` with a floor value, for example `0`,
in which case `limma` will still include this entry in the contrast.

## bug fixes

* `se_normalize()` was updated to handle assignment to `assays(se)`
while also defining rows and columns, otherwise a bug that creates
assay name `NA`; now properly assigns the assay name.

# slicejam 0.0.8.900

## new functions

* `genomic_regions_from_gtf()` - stand-alone function to convert
a GTF file into `genome_regions`, in the form of a
`GRanges` objects annotated for `"Promoters", "TTS", "exons", "introns"`,
which are also annotated by the gene defined in the GTF file.
A custom subset of `detectedTx` or `detectedGenes` can be supplied,
representing only the active or detected genes relevant to the
experiment.
* `se_normalize()` provides several normalization methods, which can
be applied to a `SummarizedExperiment` object. Multiple normalizations
can be applied, each is stored in the output `SummarizedExperiment`
object.
* `matrix_normalize()` is the base function for normalizing a numeric
data matrix. It operates on one numeric matrix, and performs
only one normalization method, from a selection of available methods.
* `update_list_elements()` and `update_function_params()` are two helper
functions, used to combine default parameters from a list, with a subset
of customized parameters. It is used to send normalization parameters
in a list named by one or more normalization methods.
* `statsdf2bed()` is a slicejam-specific function, used to convert the
stats `data.frame` into BED format for export to a file.

# slicejam 0.0.7.900

## enhancements

* `slicejam_analysis.Rmd` output indicates when sections are skipped:

    * Venn diagrams are skipped when there is only one contrast.
    * Heatmaps are not created when there are no MGM hits.
    * ATAC-mode is turned off when `ATAC=0`.
    * RData file is not saved when `SAVE_RDATA=0`.

* `slicejam_analysis.Rmd` has new argument `CURATION_TXT` which
is a tab-delimited text format. This mechanism has
replaced `CURATION_YAML` as more scientist-friendly. Example
format:

         Pattern            Batch         Group
         NOV14_p2w5_VEH     NOV14         p2w5_Veh
         NOV14_p4w4_VEH     NOV14         p4w4_Veh
         NOV14_UL3_VEH      NOV14         UL3_Veh
         NS644_UL3VEH       NS644         UL3_Veh
         NS50644_UL3VEH     NS50644       UL3_Veh
         NS644_p2w5VEH      NS644         p2w5_Veh

* Removed `CURATION_YAML` argument altogether.
* `slicejam_analysis.Rmd` now always exports two files:

    1. `"curated_samples.txt"` which contains the table of sample annotations
    and filenames used in the analysis.
    2. `"contrasts.txt"` which contains each statistical contrast one
    per row. If this file is edited, the edited contrast names will be
    used for the analysis.
    
* Suggested above, `"contrasts.txt"` can be supplied with a specific list
of statistical contrasts, which will be used in downstream analysis.
Currently the default `maxDepth=1` limits contrasts to one-way comparisons,
however the function `limma::makeContrasts()` is used to recognize
contrast names, which means you can encode your own two-way contrasts,
using format: `(group_one-group_two)-(group_three-group_four)`. Make
sure the group names match the values in the file `"curated_samples.txt"`.

## New R functions

* `curate_to_df_by_pattern()` is a function that handles the `curation.txt`
parsing. The function name may change in future.
* `import_featurecounts()` is a simple import function for featureCounts
text output files. In future this function may also apply curation by
`curation.txt` logic.

# slicejam 0.0.6.900

## enhancements

* `how_to_slicejam.md` was updated with much more comprehensive
documentation of the workflow, including how to create `"curation.yaml"`,
and quick start guide to run the pipeline.
* `TODO.md` was organized with more current set of ideas:
annotate BroadPeaks to help filter for multi-slice hits from
the same BroadPeak; user-defined "Problem_Regions.bed" to filter
problem peak slices before normalization (e.g. centromeres,
repeats, etc.)

# slicejam 0.0.5.900

## enhancements

* `slicejam_analysis.Rmd` updated to fix numerous small text formatting
issues; fixed output of inline R code (`r cat(X)` does not work,
but `r X` does work, go figure). Fixed incidence matrix kable output
(`print(kable_coloring(...))` does not work, but `kable_coloring(...)`
does.

## changes to existing analysis steps

* The MA-plots created by `jamma::jammaplot()` now use x-axis
minimum `1` instead of previous x-axis `4` -- the previous value
reflected typical lowest max group mean (MGM) value `4`, and
helped visually focus on data above this threshold. The new baseline
avoids noise near zero (x-axis `2` in log2(1+x) space refers to 3
read counts). Also, a vertical line is drawn on each plot showing
the current MGM threshold.

## enhancements

* `slicejam_analysis.Rmd` does better at checking for `curation_yaml`
and properly falling back to manual parsing when the yaml step fails.
* Colors for the design and contrast matrix tables were changed to
`-1="blue"` and `1="red"`, they got switched at some point.

# slicejam 0.0.3.900

## enhancements

* BED format actually matches BED format (doh)
* BED files include a track header line, but since the track name
can only be 15 characters, 12 character are taken by year, month, day,
hour, minute.
* Renamed output bed from from "statshits" to "stathits" avoiding
embarrassing unintended substrings.

# slicejam 0.0.2.900

## enhancements

* `GROUPCHECK=1` environment variable will optionally export
only the sample group information and stop before performing
any data analysis.
* BED files are exported: all peaks; MGM-filtered peaks; and MGM-filtered
differential peaks.
* Added `feature_type_winner` to the BED output name.
* Created `"how_to_slicejam.md"` to describe how to run the workflow
using Rscript along with `"run_slicejam.R"`.

# slicejam 0.0.1.900

## bug fixes

* Modified `slicejam_analysis.Rmd` section
`"Peak to ATAC Promoter Analysis"`
to use `ignore.strand=TRUE` while overlapping ATAC peaks by
`genome_regions`.
* Modified `slicejam_analysis.Rmd` to use `select="all"` when
annotating peaks to nearest gene, which allows peaks to be
annotated to one or more genes that have the same lowest
distance.

## changes to slicejam_analysis.Rmd

* BED files are created for three subsets of peaks:

   * all peaks
   * MGM-filtered peaks
   * MGM-filtered peak hits

# slicejam 0.0.0.900

* Initial package creation.
