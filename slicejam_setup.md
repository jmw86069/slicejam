# How to set up slicejam

This document describes steps required to set up slicejam
for analysis.

## Dependencies

This Rmarkdown currently requires several R packages which
are included in the package DESCRIPTION file.

In addition, this workflow requires functions from an
external R library of custom functions `"CBioRUtils.R"`
which is not yet included in this or any other R package.
(Work in progress.) Contact James M. Ward for access
to this file, if necessary.

## Rmarkdown

The Rmarkdown file `slicejam_analysis.Rmd` is intended to be run
as a command line tool. As such, it has several configuration
options to guide the workflow.

### Options via Environment Variables

Environment variables are used to pass most configuration options
to the Rmarkdown script.

* `OUTDIR` - output directory, by default the local directory `"."`
* `CACHEDIR` - the directory to store cache files. By default
it creates a subfolder under the `OUTDIR` directory, for example
`"./cache"`. The cache folder contains all the images used in
the summary document, in PNG and PDF format, so they can be used
as figures in other presentations.
* `SAVE_RDATA` - logical (1|0) indicating whether to store the RData file
with the full R session. Use `SAVE_RDATA=1` to save RData, or
anything else like `SAVE_RDATA=0` not to save RData. This RData
file can be pretty large, but is also sometimes useful to
review intermediate data.
* `ATAC` - logical (1|0) indicating whether the input data contains
ATAC-seq data. When ATAC-seq data is used, the peaks are annotated
differently than with general ChIP-seq data.
* `GTF` - full file path or URL to a GTF file that will be used to
annotate peaks. This GTF file is also used to define genome regions
such as "promoter", "exon", "intron", "TTS", and "intergenic".
* `FC_FILE` - the featureCounts output file, usually with `".fc"`
file extension.
* `CURATION_YAML` - a full file path or URL to a `yaml` file
used for curating the featureCounts column header into
statistical grouping. See below for more details.
* `MGM` - max group mean threshold used to filter peaks included
in the differential peak statistical testing. This filter requires
the mean peak signal to be at least this threshold in at least one
sample group. If no group has signal that meets this threshold,
it is discarded from the statistical analysis as if it were a
null region. It vastly improves the adjusted P-value, by reducing
the number of peaks tested.


### Other Helpful Environment Variables

* `R_LIBS_USER` - the path or paths to R package libraries. Using
a specific path lets one or more users share the same R package
installation directory, which is helpful if you ever run the Rmarkdown
as a different user than the user that installed the required R
packages.

### Statistical Design and Sample Groups

A critical step during the analysis is to define appropriate
statistical design, in this case in terms of sample groups.
There are three key pieces of information used by this workflow:

1. Group: Each column in the featureCounts file should be assigned
to one group. A group can have one or more members, and should represent
something close to a biological or technical replicate of the
same experimental perturbation.
2. Run: A "Run" is defined as a sequencing run, or batch, and in
fact when there are multiple "Run" values in an experiment,
batch adjustment is performed to correct any effects between the
two runs. In many cases there may only be one "Run".
3. Sample: Effectively an abbreviated label that uniquely
identifies each column in the featureCounts file.

Ideally, the Group, Run, and Sample can be derived from the column
header in the featureCounts file. But reality is what it is, and
sometimes the column header is not perfectly formatted. For this
purpose we have a curation workflow.

#### Curation of featureCounts Columns to Design Elements

The most configurable method to convert featureCounts column headers
to group information is to use a `"curation.yaml"` file. The
`"splicejam"` R package has a function `curateVtoDF()` which
takes a vector (in this case column headers from featureCounts),
and applies curation rules to produce a `data.frame` output.

An example `"curation.yaml"` file:

```
Group:
- - WT|wildtype
  - WT
- - Mut|mutant
  - Mut
Run:
- - NS[0-9]+
  - \\1
Sample:
- - [^/]+$
  - \\1
```

The basic curation workflow:

* It creates a column `"Group"`, where any entry
that matches the regular expression `"WT|wildtype"`
(matches "WT" or "wildtype") is assigned the value `"WT"`.
Any entry that matches `"Mut"` or `"mutant"` is
assigned `"Mut"`. Everything else is assigned the original value,
which should usually be avoided.
* It then creates a column `"Run"` where any entry that matches
`"NS[0-9]+"` will be assigned that matched value. For example
`"NS50077.SampleA"` will match `"NS50077"`.
* It then creates a column `"Sample"` which matches
`"[^/]+$`, which matches any character except `"/"`.
For example `"/some/file/folder/NS50077.SampleA.bam"` will
match `"NS5077.SampleA.bam"`.

Note that values in the `"Sample"` column are forced to be
unique by calling `jamba::makeNames()`. When values are already
unique, they are not changed, but any duplicated values
are given a suffix `"_v1"`, `"_v2"`, etc.

The `"Sample"` column is used as the label on several output
plots such as the MA-plots. The `"Group"` label is used
also, but not as the per-replicate label.

#### Fallback sample groups

In the absence of `"CURATION_YAML"` environment variable,
or if the file does not exist, the format is assumed to
be something like this:

**RUN**_*GROUP*_**REP**

* Samples are defined by calling `basename()` to remove the
directory of each file (if it exists), then removing
the file stem, and any recognized keywords such as `"dedup"`,
`"cutadapt"`, `"fastq"`, etc.
* Sample is split using the delimiter `"[-_.]"` and the
first word is used for `"Run"`, the second word is used
for `"Group"`.
* The delimiter therefore can be any of: hyphen `"-"`,
underscore `"_"`, or decimal `"."`.
* Group is converted to factor, whose levels are ordered
so that several recognized control labels are first:
`"Control","ctrl","ctl","veh","blank","ntc"`, followed by
the groups in the order they appear as column names in
the featureCounts file.

### GTF processing

The GTF file provided is processed before used during
peak annotation. After processing, intermediate files
are stored (if possible) for faster re-use.
The file is saved with extension `.genome_regions.RData`.

If the `.genome_regions.RData` cannot be found in the
same directory as the GTF, it will look in the current
working directory (where Rmarkdown is invoked).

Similarly, it the output file cannot be written to
the same directory as the GTF file, it will attempt to
write in the current working directory (where Rmarkdown
is invoked).
