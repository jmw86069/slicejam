usage

export GTF=some_gtf_file.gtf
RHOME=/path/to/R_HOME

set RHOME with
export RHOME=$(R361 --vanilla --slave -e 'Sys.getenv("R_HOME")')

export R_LIBS_USER=/ddn/gs1/home/wardjm/R/x86_64-pc-linux-gnu-library/3.6


## then to run the script for a BED file of peaks:
/ddn/gs1/home/jmw86069/Projects/server_scripts/genomic_regions_from_gtf --gtf=${GTF} --bed=some_bed_file.bed
