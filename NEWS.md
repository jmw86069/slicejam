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
