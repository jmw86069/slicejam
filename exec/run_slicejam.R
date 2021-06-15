## Rscript
##
## run_slicejam.R
##
## small wrapper around gokeyATAC_analysis.Rmd
##
## Run like this:
##
## GTF=genome.gtf ATAC=1 FC_FILE=071419_2124_38to100_NoY_NoM_NoProblems.fc.txt CURATION_YAML=curation.yaml DRYRUN=1 Rscript run_gokeyATAC.R
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
## CURATION_YAML=curation.yaml optional YAML file
##    to define columns: Group, Sample, Run using the featureCounts colnames
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

sliceargs <- get_slicejam_args();

DRYRUN <- sliceargs$DRYRUN;
ATAC <- sliceargs$ATAC;
fc_file <- sliceargs$fc_file;
fc_filepath <- sliceargs$fc_filepath;
CURATION_TXT <- sliceargs$CURATION_TXT;
GTF <- sliceargs$GTF;
GTFNAME <- sliceargs$GTFNAME;
DRYRUN <- sliceargs$DRYRUN;
OUTDIR <- sliceargs$OUTDIR;
ATAC <- sliceargs$ATAC;
MGM <- sliceargs$MGM;
NORM <- sliceargs$NORM;
NORMMIN <- sliceargs$NORMMIN;
NORMSHORT <- sliceargs$NORMSHORT;
fc_base <- sliceargs$fc_base;
fc_basedir <- sliceargs$fc_basedir;
fc_html <- sliceargs$fc_html;
knit_root_dir <- sliceargs$knit_root_dir;
cache_dir <- sliceargs$cache_dir;

if (FALSE) {
   DRYRUN <- Sys.getenv("DRYRUN");
   fc_file <- Sys.getenv("FC_FILE");
   CURATION_TXT <- Sys.getenv("CURATION_TXT");
   MGM <- Sys.getenv("MGM");
   if (length(MGM) == 0 || nchar(MGM) == 0) {
      MGM <- 4;
   }
   
   OUTDIR <- Sys.getenv("OUTDIR");
   
   GTF <- Sys.getenv("GTF");
   if (length(GTF) > 0 && nchar(GTF) > 0) {
      if (!file.exists(GTF)) {
         stop(paste0("GTF file not found:", GTF));
      }
      GTFNAME <- Sys.getenv("GTFNAME");
      if (length(GTFNAME) == 0 || nchar(GTFNAME) == 0) {
         GTFNAME <- paste0("_",
            gsub("gene[s]*", "",
               gsub("[-_ ]+", "",
                  gsub("[.][^.]+$", "", basename(GTF)))));
      }
      gtf_stem <- paste0("_gtf", GTFNAME);
      Sys.setenv(GTFNAME=GTFNAME);
   } else {
      stop("GTF must be supplied.");
      gtf_stem <- "";
   }
   
   ATAC <- Sys.getenv("ATAC");
   if (length(ATAC) == 0 || nchar(ATAC) == 0 || ATAC %in% c("0")) {
      ATAC <- 0;
   } else {
      ATAC <- 1;
   }
   
   NORM <- Sys.getenv("NORM");
   if (length(NORM) == 0 || nchar(NORM) == 0 || grepl("quant", ignore.case=TRUE, NORM)) {
      NORM <- "quantile";
      NORMSHORT <- "quant";
   } else if (grepl("med.*gr.*p", ignore.case=TRUE, NORM)) {
      NORM <- "mediangroup";
      NORMSHORT <- "medgrp";
   } else if (grepl("med.*", ignore.case=TRUE, NORM)) {
      NORM <- "median";
      NORMSHORT <- "med";
   } else if (grepl("none", ignore.case=TRUE, NORM)) {
      NORM <- "none";
      NORMSHORT <- "none";
   } else {
      NORM <- "quantile";
      NORMSHORT <- "quant";
   }
   
   if (nchar(fc_file) == 0) {
      stop("FC_FILE is not defined.");
   }
   if (!file.exists(fc_file)) {
      stop(paste0("FC_FILE is not found:", fc_file));
   }
   
   fc_base <- gsub("[.]fc$|[.]fc[.]txt$", "",
      fc_file);
   ## add the max group mean threshold (mgm) _mgm4
   fc_base <- paste0(fc_base, "_mgm", MGM);
   
   ## add the ATAC mode, _atac1 when ATAC-mode is enabled
   if (!ATAC %in% c(0)) {
      fc_base <- paste0(fc_base, "_atac", ATAC);
   }
   
   ## Add GTF stem using GTFNAME or small name _gtfhg19121619
   if (nchar(OUTDIR) == 0) {
      fc_base <- paste0(fc_base, gtf_stem);
      fc_basedir <- paste0(fc_base,  "_analysis");
   } else {
      fc_base <- OUTDIR;
      fc_basedir <- OUTDIR;
   }
   
   fc_html <- paste0(fc_base,  ".html");
   
   knit_root_dir <- file.path(
      getwd(),
      fc_basedir);
   fc_filepath <- normalizePath(fc_file);
   
   cache_dir <- file.path(fc_basedir, "cache", "");
}

jamba::printDebug("FC_FILE:      ", fc_file);
jamba::printDebug("output_dir:   ", fc_basedir);
jamba::printDebug("output_file:  ", fc_html);
jamba::printDebug("knit_root_dir:", knit_root_dir);
jamba::printDebug("cache_dir:    ", cache_dir);
jamba::printDebug("fc_file path: ", fc_filepath);
jamba::printDebug("CURATION_TXT: ", CURATION_TXT);
jamba::printDebug("GTF:          ", GTF);
jamba::printDebug("GTFNAME:      ", GTFNAME);
jamba::printDebug("ATAC:         ", ATAC);
Sys.setenv(FC_FILE=fc_filepath);
Sys.setenv(CACHEDIR=cache_dir);

if (!"0" %in% DRYRUN) {
   jamba::printDebug("DRYRUN mode.");
   stop("Exiting due to DRYRUN=1, set DRYRUN=0 to proceed.");
}

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
