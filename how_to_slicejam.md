# How to slicejam

## Quickstart guide to run slicejam:

The examples here are documented in detail in other sections
of this guide.

Before running slicejam you need these elements:

* featureCounts .fc file
* GTF file
* curation.yaml file
* Path to `Rscript` for the correct R version

The sections describe how to create these files.

### Check that the input data is correct

#### 1. Run slicejam for non-ATAC data, with `DRYRUN=1` mode to check that all files exist.

```
DRYRUN=1 \
   GTF=hg19_gencode.gtf \
   FC_FILE=your_favorite_featurecounts_file.fc \
   CURATION_YAML=curation.yaml \
   Rscript run_slicejam.R
```

#### 2. Run slicejam with `GROUPCHECK=1` to verify the `curation.yaml` defines correct sample groups.

```
DRYRUN=0 \
   GROUPCHECK=1 \
   GTF=hg19_gencode.gtf \
   FC_FILE=your_favorite_featurecounts_file.fc \
   CURATION_YAML=curation.yaml \
   Rscript run_slicejam.R
```

### Run the full pipeline

#### 3. Run slicejam for non-ATAC data, with `DRYRUN=0`:

```
DRYRUN=0 \
   GTF=hg19_gencode.gtf \
   FC_FILE=your_favorite_featurecounts_file.fc \
   CURATION_YAML=curation.yaml \
   Rscript run_slicejam.R
```

#### 4. Run slicejam with `ATAC=1` and `DRYRUN=0` to run the pipeline.

```
DRYRUN=0 \
   GTF=hg19_gencode.gtf \
   ATAC=1 \
   FC_FILE=your_favorite_featurecounts_file.fc \
   CURATION_YAML=curation.yaml \
   Rscript run_slicejam.R
```

## Prerequisites

1. featureCounts output file with `".fc"` file extension. (As created by `mergeSplitCountBed`)

   * Any `.fc` file will work, but pay attention to column headers
   for each sample.
   * The column headers will be used to define sample groups in Step 2.

2. Create a `"curation.yaml"` file (see below).

   * The simplest form converts a sample ID to `"Group"` (to define sample
   replicates per group), `"Run"` (sequencing run, used for batch adjustment),
   and `"Sample"` (used as a short visual label).

3. GTF file that defines gene-transcript-exons.

   * Column 9 is expected to contain at least these attributes:
   
      * `"gene_name"`
      * `"gene_id"`
      * `"transcript_id"`
      
   * For example, human hg19 gencode v32 is located here:
   
         /ddn/gs1/shared/dirib/reference_genomes/hg19/hg19gencode/gencode.v32lift37.annotation/gencode.v32lift37.annotation.gtf
         
   * You can copy the GTF file to your analysis directory, or
   you can point to that file.
   * You can be fancy and define a variable to store the file path
   
         export GENCODE_GTF_V32="/ddn/gs1/shared/dirib/reference_genomes/hg19/hg19gencode/gencode.v32lift37.annotation/gencode.v32lift37.annotation.gtf"

4. Two files from the slicejam R package.

   * `"run_slicejam.R"`
   * `"slicejam_analysis.Rmd"`

5. Path to the correct `Rscript` for the correct version of R, we recommend R-3.6 or higher.

   * For NIEHS servers, `"R361"` is used for R-3.6.1. The path
   can be found by running this command:
   
         ls -l $(which R361)
            
      * The output shows:
      
            /ddn/gs1/biotools/misc/bin/R361 -> /ddn/gs1/biotools/R361/bin/R

         which means the actual `R` is located at
      
            /ddn/gs1/biotools/R361/bin/R

      * Edit this path to change R to Rscript:
      
            /ddn/gs1/biotools/R361/bin/Rscript
            
      * Confirm the file exists by running
      
            ls -l /ddn/gs1/biotools/R361/bin/Rscript
            
   * The same steps can be run for R:

         ls -l $(which R)
         
   * Confirm the version of R is 3.6 or higher with the argument `"--version"`:
   
         /ddn/gs1/biotools/R361/bin/Rscript --version

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

## slicejam options

Slicejam has the following options:

   * `DRYRUN`
      - use `DRYRUN=1` to run dry-run mode
      - use `DRYRUN=0` to run the full analysis
      
   * `FC_FILE` - the featureCounts `.fc` input file
      - for example: `FC_FILE=some_peak_counts.fc`
      
   * `GTF` - file path to the GTF file
      - for example: `GTF=gencode.v32lift37.annotation.gtf`
      
   * `GTFNAME` - optional label
      - for example: `GTFNAME=gencodev32`
   
   * `MGM` - max group mean value
      - for example: `MGM=4` applies a cutoff to normalized
         group mean counts which were transformed by `log2(1+x)`.
         The cutoff therefore requires `15` or
         more normalized counts in at least one sample group.
      - for example: `MGM=6` may be best for broader peaks,
         1kb or larger.
      - The `MGM` value can be informed by the MA-plots,
         after running slicejam once with default values.
      
   * `CURATION_YAML`
      - path to a `"curation.yaml"` file, used to convert
      featureCounts column names to sample groups. See below.
      - for example: `CURATION_YAML=curation.yaml`
      - if no `CURATION_YAML` is defined, all samples will
         be processed as individual samples without grouping.

   * `GROUPCHECK`
      - use `GROUPCHECK=1` so the pipeline will only print
      sample groups defined by `curation.yaml`. This step is
      good to run the first time, to verify correct sample
      grouping.
      - use `GROUPCHECK=0` to run the full pipeline.
      
   * `ATAC`
      - use `ATAC=1` to enable ATAC-mode, which annotates
      peaks by only promoters that overlap peaks above the `MGM`
      threshold.
      - use `ATAC=0` to disable ATAC-mode.
      - In both cases, peaks are also annotated by nearest
      gene body defined in the `GTF` file. ATAC-mode additionally
      annotates peaks by nearest TSS of genes whose promoters
      overlap a MGM-filtered peak.

### Options can be defined one of two ways.

1. Define options in one command.

```
DRYRUN=1 FC_FILE=feature_counts.fc GTF=gencode_hg19.gtf CURATION_YAML=curation.yaml Rscript run_slicejam.R
```

This command can be split across multiple lines by adding `"/"`
to the end of each line, to make it more readable:

```
DRYRUN=1 \
   FC_FILE=feature_counts.fc \
   GTF=gencode_hg19.gtf \
   CURATION_YAML=curation.yaml \
   Rscript run_slicejam.R
```

2. Define variables before running slicejam.

The variables are each defined using `"export"` as shown:

```
export DRYRUN=1
export FC_FILE=feature_counts.fc
export GTF=gencode_hg19.gtf
export CURATION_YAML=curation.yaml
```

You can confirm the value of a variable with the command `"echo"`:

```
echo ${GTF}
```

Next run slicejam:

```
Rscript run_slicejam.R
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


## How to interpret the HTML report

### Review peak distribution

The first plot shows the distribution of peak/feature widths
using genome coordinates.

### Confirm correct sample group design

Samples are assigned to sample groups. Confirm that this step
resulted in the intended experiment design.

### Review MA-plots

MA-plots are used to review data quality, and data normalization.
Most downstream statistical methods make assumptions of
the data, and many of those assumptions can be visualized
in MA-plots.

#### MA-plot progression

MA-plots are produced for progressively specific scenarios, starting
with the data overall.

* **Global-centered raw data**: The data is directly compared across
all sample replicates. If there are batch effects, or substantial
differences in read depth, you should see it here.
* **Global-centered quantile-normalized data**: The data is
quantile-normalized then compared across all samples. If there are
no batch effects, the data should look relatively clean at this step.
* **Group-centered quantile-normalized data**: The normalized data
is centered within each sample group, which means the MA-plots
should represent only the difference from group average for each
sample replicate. Ideally, the data should be centered at y=0
with low variance.
* **Group-centered, batch-adjusted, quantile-normalized**: The
normalized data is batch-adjusted across sequencing runs, 
which only occurs when there are multiple `"Run"` values.
If there is a batch effect, the variation here should be substantially
lower than the previous MA-plot which is not batch-adjusted.

#### MA-plot for sample outliers

MA-plot panels are colored yellow when its median variance is
more than two-tiemes higher than the median variance across
all samples. This threshold is a 2x MAD factor cutoff, which
roughly means the variance is more than twice the variance
of other samples. This typically represents a technical outlier,
and not a biological outlier, and therefore would adversely
affect the statistical modeling of true biological variance
in the experiment.

We have chosen not to automate the process of sample outlier
detection and removal, and leave this decision-making to
scientists. The key decision is whether a sample represents a
technical failure, in which case including its data would
not reflect the true underlying biology.

To remove an outlier sample, re-run
slicejam with a `.fc` file that does not include the outlier sample.

#### MA-plot for statistical thresholds

MA-plots represent variance as y-axis scatter, compared
to the average signal on the x-axis. The variance is usually
much higher with low signal, since the shot noise (imposed by
getting integer counts for low number of counts) becomes higher
than biological/technical noise. The efficiency of the IP can
also impose variance in the form of high background signal,
which produces a wide variance at low counts, dependent
upon the peak width.

Upon reviewing the MA-plots, it may become clear what threshold
should be used to be above this "shot noise" or "background noise",
which helps define a usable subset of "detected peaks."

We have avoided automating this baseline detection, for
now the threshold is user-defined, with default `MGM=4`
translating to 15 or more read counts per peak. For wide
peaks or higher background signal, `MGM=6` (31+ counts)
or `MGM=7` (63+ counts) are suitable.


### Review sample correlation heatmaps

Samples are centered per row, which subtracts the average
signal for each peak, producing data that contains the difference
from average for each peak. In general, there should be no correlation
if the remaining data has more noise than signal. However, when a
sample differs from average in the same way as another sample, they
share positive correlation.

In general we expect sample to have higher correlation to samples
of the same sample group, and lower correlation to other sample
groups.

#### Sample correlation of quantile-normalized data

* If the data has a batch effect, the heatmap typically clusters
batch together, then sample group.
* If there is no batch effect, samples should typically cluster
together based upon sample group.

#### Sample correlation of quantile-normalized, batch-adjusted data

* If the data has a batch effect, this heatmap should cluster
samples together based upon sample group, and not by sequencing run.
In other words, this heatmap should better represent the experiment
design than the previous correlation heatmap.
* If the data has no batch effect, this heatmap should closely
resemble the previous correlation heatmap.

> In general, applying batch adjustment to data with no batch
effect, has very little impact on the data. That impact is measureable,
so if needed, the data can be directly compared.


