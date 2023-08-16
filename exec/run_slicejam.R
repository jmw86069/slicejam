## Rscript
##
## run_slicejam.R
##
## small wrapper around gokeyATAC_analysis.Rmd
##
## Run like this:
##
## GTF=genome.gtf ATAC=1 FC_FILE=071419_2124_38to100_NoY_NoM_NoProblems.fc.txt CURATION_TXT=curation.txt DRYRUN=1 Rscript run_gokeyATAC.R
##
## GTF="genome.gtf" # the GTF file with gene-exon features
##
## GTFNAME="genome_label" # optional label used in place of GTF file stem
##
## MGM defines an optional maxGroupMean threshold, default MGM=4.
##
## ATAC=1 defines ATAC mode which calls active promoters that directly
##    overlap ATAC peaks, thus reducing the genes used during "nearest gene"
##
## DRYRUN=1 will run in dry-run mode, printing only what it will do,
##    but will not run the Rmarkdown analysis.
##
## FC_FILE=featureCounts_output.fc
##
## CURATION_TXT=curation.txt optional tab-delimited text file
##    with columns: Filename, Sample, Run, Group, Label (optional)
##    where Filename should match some substring of the column headers
##    in the featureCounts file
##
## It will remove the '.fc' or '.fc.txt' file extension, and
## create a file with extension '_analysis.html' with output
## from the Rmarkdown workflow.
##
## It will place figures in a subfolder with the
## base filename plus the extension '_analysis'
##
## Figures will include PNG images rendered in the HTML file, and
## PDF versions of each figure, for more accurate editing.

sliceargs <- slicejam::get_slicejam_args(verbose=TRUE);

DRYRUN <- sliceargs$DRYRUN;
ATAC <- sliceargs$ATAC;
fc_file <- sliceargs$fc_file;
fc_filepath <- sliceargs$fc_filepath;
CURATION_TXT <- sliceargs$CURATION_TXT;
GTF <- sliceargs$GTF;
GTFNAME <- sliceargs$GTFNAME;
OUTDIR <- sliceargs$OUTDIR;
MASK <- sliceargs$MASK;
MGM <- sliceargs$MGM;
NORM <- sliceargs$NORM;
NORMMIN <- sliceargs$NORMMIN;
NORMSHORT <- sliceargs$NORMSHORT;
fc_base <- sliceargs$fc_base;
fc_basedir <- sliceargs$fc_basedir;
fc_html <- sliceargs$fc_html;
knit_root_dir <- sliceargs$knit_root_dir;
CACHEDIR <- sliceargs$CACHEDIR;

if (DRYRUN) {
   jamba::printDebug("DRYRUN mode.");
   stop("Exiting due to DRYRUN=1, set DRYRUN=0 to proceed.");
}

# Define environment variables updated with full path
Sys.setenv(FC_FILE=fc_filepath);
Sys.setenv(OUTDIR=OUTDIR);
Sys.setenv(CACHEDIR=CACHEDIR);
Sys.setenv(CURATION_TXT=CURATION_TXT);
Sys.setenv(GTF=GTF);

# confirm slicejam_analysis.Rmd is in the current directory
Rmd_file <- system.file(package="slicejam", "exec/slicejam_analysis.Rmd");
if (!file.exists("slicejam_analysis.Rmd")) {
   jamba::printDebug(c("Copying ",
      "slicejam_analysis.Rmd",
      " from slicejam package into current directory."),
      sep="");
   file.copy(from=Rmd_file,
      to="slicejam_analysis.Rmd");
}

if (!file.exists("slicejam_analysis.Rmd")) {
   stop("slicejam_analysis.Rmd file does not exist in the current directory.");
}
rmarkdown::render("slicejam_analysis.Rmd",
   output_format="html_document",
   knit_root_dir=knit_root_dir,
   output_dir=knit_root_dir,
   intermediates_dir=knit_root_dir,
   output_file=fc_html
)
