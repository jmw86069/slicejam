# Steps to use slicejam

## Prerequisites

1. featureCounts output file with `.fc` file extension.
2. For sample grouping, create a `"curation.yaml"` file.
3. Provide a GTF file that contains transcript exons.

   * Column 9 is expected to contain `"gene_name"`, `"gene_id"`,
   and `"transcript_id"`.

4. Copy two files from the slicejam R package to your analysis directory:

   * `"run_slicejam.R"`
   * `"slicejam_analysis.Rmd"`

5. Find `Rscript` appropriate for your version of R.

   * Some linux hosts have multiple versions of R installed,
   find the path to the correct one if this is the case.
   * For example, in our linux environment I use `R361` to run R version 3.6.1:
   
      * `ls -l $(which R361)`
      * For our environment, the result shows that `R361` is a symbolic link:
      
      > `/ddn/gs1/biotools/misc/bin/R361 -> /ddn/gs1/biotools/R361/bin/R`
      
      * Follow the symbolic link to find `Rscript` for this version of R.
      * `ls -l /ddn/gs1/biotools/R361/bin/R*`

## Run slicejam using Rscript

1. The `"run_slicejam.R"` file is a small wrapper script that calls
`rmarkdown::render()` in the proper way.
2. Define environment variables for your workflow:

   * `FC_FILE` - the featureCounts `.fc` input file
   * `DRYRUN` - use `DRYRUN=1` to run dry-run mode, `DRYRUN=0` to run the full analysis
   * `GTF` - file path to the GTF file
   * `GTFNAME` - optional label to use for the GTF file
   * `MGM` - max group mean value, usually `MGM=4` or higher, where `4`
   refers to the cutoff for data transformed by `log2(1+x)`, and therefore
   means 15 or more normalized counts in at least one sample group mean.
   * `CURATION_YAML` - path to a `"curation.yaml"` file, to convert
   featureCounts column names to sample groups.
   * `ATAC` - if your data contains ATAC-seq peaks use `ATAC=1`; otherwise `ATAC=0`

3. Call `Rscript run_slicejam.R`.

   * If using a linux/Mac machine, make sure to `export` environment variables,
   for example `export DRYRUN=1`.

## By-products of Rmarkdown output

1. The Rmarkdown file is used to produce an HTML summary report of the analysis.
2. The HTML output does not contain the images, but points to a sub-folder.
The intent was to allow re-using the images directly, also to save a PNG and
PDF format at the same time. The PDF format preserves detail even when resized,
the PNG format is composed of image pixels that do not scale well.
3. A subfolder is created which contains text `.txt` stats files, and BED `.bed`
peak files.
4. A subfolder is created which contains a cache folder, used for Rmarkdown
caching in a specific way to prevent sharing this cache with other analysis
runs.

   * Typically, Rmarkdown creates a cache subfolder based upon the source
   Rmarkdown file, therefore all analyses that use the same Rmarkdown file
   also share the same cache folder -- which is not good for this workflow.
   In that scenario, no cache files will ever be re-used.


## How to create a curation.yaml file

An important and useful part of the workflow is converting column
names in the featureCounts `.fc` file into useful sample grouping.
The `"curation.yaml"` file defines rules to create three values
for each replicate:

1. Sample - some label or unique identifier for each sample replicate
2. Group - group name for each sample replicate
3. Run - the sequencing run, used when batch-adjustment should be performed

The `"curation.yaml"` file should have three sections, one
each for `"Sample:"`, `"Run:"`, and `"Group:"`. The content is
described in the context of an example below:

      Run:
      - - ^[A-Z]+[0-9]+
        - \\1
      Group:
      - - (p[0-9]w[0-9]|UL3)
        - \\2
      Sample:
      - - ((p[0-9]w[0-9]|UL3).*_v[0-9]+)
        - \\2

For the section `Run:` there is one rule:

* It matches `"^[A-Z]+[0-9]+"` - which matches uppercase
letters, then numbers, starting at the beginning of the
column name. If the pattern is matched, that pattern is
used as the value.
* Therefore, `"NOV215_p2w5_blahblah.bam"` is converted to `"NOV215"`

For the section `Group:` there is one rule:

* It matches `"(p[0-9]w[0-9]|UL3)"` which matches either
`"p[0-9]w[0-9]"` or `"UL3"`.
* Therefore, `"NOV215_p2w5_blahblah.bam"` is converted to `"p2w5"`.
* Similarly, `"NOV215_UL3_blahblah.bam"` is converted to `"UL3"`.

For the section `Sample:` there is one rule:

* It matches `"((p[0-9]w[0-9]|UL3).*_v[0-9]+)"`. This rule extends
the previous rule `"(p[0-9]w[0-9]|UL3)"` which matches either
`"p[0-9]w[0-9]"` or `"UL3"`. The rule adds `".*_v[0-9]+"` which
matches any string until it reaches `"_v"` followed by one or more
numbers.
* Therefore, `"NOV215_p2w5_blahblah_v1.bam"` is converted to `"p2w5_blahblah_v1"`.
* Similarly, `"NOV215_UL3_blahblah_v2.bam"` is converted to `"UL3_blahblah_v2"`.

There can be multiple rules per column if necessary. Rules can be simpler,
for example.

      Group:
      - - p2w5
        - p2w5
      - - p4w4
        - p4w4
      - - UL3
        - UL3

There are three rules, and in each case the pattern match is also used
as the replacement. This technique is effective for curating inconsistencies,
for example:

      Group:
      - - p2w5DEX
        - p2w5_dex
      - - p2w5_DEX
        - p2w5_dex

