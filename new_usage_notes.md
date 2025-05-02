# usage of genomic_regions_from_gtf

export GTF=some_gtf_file.gtf
RHOME=/path/to/R_HOME

set RHOME with
`export RHOME=$(R361 --vanilla --slave -e 'Sys.getenv("R_HOME")')`

Optionally you may need to point to a specific directory
in which the required R packages are installed.
`export R_LIBS_USER=/path_to/R/x86_64-pc-linux-gnu-library/3.6`

Add the path to genomic_regions_from_gtf to PATH
RPATH=/path_to/Projects/server_scripts
PATH=${RPATH}:${PATH}

Verify genomic_regions_from_gtf is accessible:
`genomic_regions_from_gtf -h`

Then run the tool on a BED file of peaks:
`genomic_regions_from_gtf --gtf=${GTF} --bed=some_bed_file.bed`

## Optional mask regions

If you have a mash file, which is a BED format file that
contains genome regions that should be removed from the
annotation, define MASK and use argument `--mask`:

```
MASK="/path_to/reference_genomes/hg19/hg19-blacklist.v2.bed"
export MASK

genomic_regions_from_gtf --gtf="${GTF}" --mask="${MASK}" --bed="${BED}"
```

## Optional detected genes

If you have a subset of genes to be used from the full GTF
annotation file, use the argument `--genes` to point to this
file.

The detected genes file should contain one column, with values
that match gene symbols in the GTF file.

```
export DETGENES="detected_genes.txt"
genomic_regions_from_gtf --gtf="${GTF}" --mask="${MASK}" --bed="${BED}" --genes="${DETGENES}"
```

When using detected genes, only the entries in the GTF
file that match these gene symbols will be used to define
the genome regions. Other entries in the GTF file are
considered "not detected" and therefore are not part
of the analysis. This option is useful when the GTF
may be comprehensive and includes a large number of
entries that are expressed only in specific cell types
which may not be relevant to the current experiment.
