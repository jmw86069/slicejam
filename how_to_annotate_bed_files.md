# How to use genomic_regions_from_gtf

## Define some environmental variables

> Note: If you are using GNU `screen`, or `tmux` remember to
> define these variables for each new shell window.

`RHOME` points the `Rscript` to the correct version of R installed
on the computer. In the example below `R361` represents the
version of R that should be used.

```
RHOME=$(R361 --no-save --slave -e 'cat(Sys.getenv("R_HOME"))')
export RHOME
```

`GTF` is the full path to a GTF file used to define genome_regions.
If this GTF file has already been used, the associated cached files
will be re-used.

```
GTF="/full/path/to/gtf/gencode.v32lift37.annotation.gtf.gz"
export GTF
```

## Ensure genomic_regions_from_gtf is on the PATH

The `genomic_regions_from_gtf` script must be on your accessible `PATH`,
which can be done one these ways:

* copy the script to a folder in your `PATH`, for example copy to `~/bin`
if that is on your `PATH`. Check path by running `echo ${PATH}`
* add the folder that contains `genomic_regions_from_gtf` to your existing
`PATH`

```
GPATH="/some/path/to/use"
PATH=${GPATH}:${PATH}
export PATH
```

## Run genomic_regions_from_gtf

### Test the configuration

```
genomic_regions_from_gtf -h
```

This step also confirms the R configuration is correct, and
will fail if the R packages or R version are incorrect.


### Define genome_regions without a BED file

You can run `genomic_regions_from_gtf` without a BED file, which
will only create or confirm that the associated genome_regions
files have been created.

```
genomic_regions_from_gtf --gtf=${GTF}
```

### Annotate BED regions using genome_regions

Finally, run `genomic_regions_from_gtf` with a BED file to annotate
BED regions with genome_regions.

```
genomic_regions_from_gtf --gtf=${GTF} --bed=BRG1KD_Weakened_peaks_q005_FC1_5.bed
```

### Annotated BED regions including MASK regions


```
MASK="/ddn/gs1/shared/dirib/reference_genomes/hg19/hg19-blacklist.v2.bed"
genomic_regions_from_gtf --gtf="${GTF}" --mask="${MASK}" --bed=BRG1KD_Weakened_peaks_q005_FC1_5.bed
```

