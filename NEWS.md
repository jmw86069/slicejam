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
