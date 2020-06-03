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

   * `DRYRUN` - use `DRYRUN=1` to run dry-run mode, `DRYRUN=0` to run the full analysis
   * `FC_FILE` - the featureCounts `.fc` input file
   * `GTF` - file path to the GTF file
   * `GTFNAME` - optional label to use for the GTF file
   * `MGM` - max group mean value, usually `MGM=4` or higher, where `4`
   refers to the cutoff for data transformed by `log2(1+x)`, and therefore
   means 15 or more normalized counts in at least one sample group mean.
   * `CURATION_YAML` - path to a `"curation.yaml"` file, to convert
   featureCounts column names to sample groups.
   * `GROUPCHECK` - use `GROUPCHECK=1` to run the pipeline only to print
   out the sample groups defined by `curation.yaml`. This option lets you
   verify the sample groups before the analysis.
   * `ATAC` - if your data contains ATAC-seq peaks use `ATAC=1` for ATAC-mode;
   otherwise `ATAC=0`. ATAC-mode annotates peaks using a subset
   of TSSes with overlapping promoter ATAC MGM-filtered peaks. These
   TSSes are "active promoter TSSes" and represent TSSes observed to be
   active in the given cell type.

3. Call `Rscript run_slicejam.R`.

   * If using a linux/Mac machine, make sure to `export` environment variables,
   for example `export DRYRUN=1`.


### Examples calling `run_slicejam.R`:

> The examples below use back-slash `\` to split the command across
multiple lines. You can also run the command on one line after removing
the back-slashes.

#### 1a. Call `run_slicejam.R` with GROUPCHECK=1 to verify the `curation.yaml`:

```
DRYRUN=0 \
   GROUPCHECK=1 \
   GTF=hg19_gencode.gtf \
   FC_FILE=071419_2124_38to100_NoY_NoM_NoProblems.fc \
   CURATION_YAML=curation.yaml \
   Rscript run_gokeyATAC.R
```

#### 2a. Call `run_slicejam.R` for non-ATAC data, in DRYRUN mode to check output:

```
DRYRUN=1 \
   GTF=hg19_gencode.gtf \
   FC_FILE=071419_2124_38to100_NoY_NoM_NoProblems.fc \
   CURATION_YAML=curation.yaml \
   Rscript run_gokeyATAC.R
```

#### 2b. Call `run_slicejam.R` for non-ATAC data, with DRYRUN mode off, to run the pipeline:

```
DRYRUN=0 \
   GTF=hg19_gencode.gtf \
   FC_FILE=071419_2124_38to100_NoY_NoM_NoProblems.fc \
   CURATION_YAML=curation.yaml \
   Rscript run_gokeyATAC.R
```

#### 3a. Call `run_slicejam.R` with ATAC-mode on, in DRYRUN mode to check output:

```
DRYRUN=1 \
   GTF=hg19_gencode.gtf \
   ATAC=1 \
   FC_FILE=071419_2124_38to100_NoY_NoM_NoProblems.fc \
   CURATION_YAML=curation.yaml \
   Rscript run_gokeyATAC.R
```

#### 3b. Call `run_slicejam.R` with ATAC-mode on, with DRYRUN mode off, to run the pipeline:

```
DRYRUN=0 \
   GTF=hg19_gencode.gtf \
   ATAC=1 \
   FC_FILE=071419_2124_38to100_NoY_NoM_NoProblems.fc \
   CURATION_YAML=curation.yaml \
   Rscript run_gokeyATAC.R
```


## By-products of Rmarkdown output

1. The Rmarkdown file is used to produce an HTML summary
report of the analysis.
2. The HTML output includes links to the images, which are stored in
a separate subfolder.

   * The images are not embedded inside the HTML, so to view the HTML
   file on another computer, copy the HTML file, and copy the subfolder
   that ends with `"_files"`.
   * The images are stored separately so they can be viewed directly.
   * All images are saved in PNG and PDF format. The PDF format preserves
   details in the plot without pixelation.

3. A subfolder is created that ends with `"_analysis"`, whose name
is derived from the input parameters.

   * The format of the folder name: `Peak_File`_`mgm4`_`GTF_name`_`files`
   * This folder contains text `.txt` stats files with all statistical
   contrasts. There is one `.txt` file with `nomgm` in the name that
   contains all peaks; and one file with `mgm4` in the name that
   contains only the peaks that met the MGM threshold.
   * This folder also contains BED `.bed` files representing:
   
      * all peaks: `_stats_nomgm_` is in the BED file name.
      * all mgm peaks: `_stats_mgm4_` is in the BED file name.
      * all differential mgm peaks: `_stathits_mgm4_` in the BED file name.

4. A subfolder is created which contains a `"cache"` subfolder,
used by Rmarkdown to retain data during processing. If the analysis
is ever interrupted, it can be re-run and it will pick up where it
left off.


## How to create a `curation.yaml` file

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
for example instead of regular expressions, they can include exact
text strings:

      Group:
      - - p2w5
        - p2w5
      - - p4w4
        - p4w4
      - - UL3
        - UL3

There are three rules, and in each case the pattern match is also used
as the replacement.
Whenever it sees `"p2w5"` anywhere in the input,
it puts `"p2w5"` in the `"Group"` column.
Whenever it sees `"p4w4"`, it puts `"p4w4"` in the `"Group"` column.
Whenever it sees `"UL3"`, it puts `"UL3"` in the `"Group"` column.
Everything else is left as-is without change, and is placed into the
`"Group"` column.

This technique is effective for curating inconsistencies,
for example we can recognize `"p2w5DEX"` and `"p2w5_DEX"`:

      Group:
      - - p2w5DEX
        - p2w5_dex
      - - p2w5_DEX
        - p2w5_dex

Most importantly, use the option `GROUPCHECK=1` to run only the
steps that `curation.yaml` steps, to verify that the output is
what you intended. When `GROUPCHECK=1` the Rmarkdown will
save a file with extension `".curated_samples.txt"` that
contains a table with the featureCounts `.fc` column headers,
and the results of the `"curation.yaml"` output.


## How to interpret the HTML report



