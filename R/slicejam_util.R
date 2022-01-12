

#' Get slicejam argument values
#' 
#' Get slicejam argument values
#' 
#' This function is a convenient wrapper to handle conversion
#' of argument values sent via environmental variables, into
#' values to be used during slicejam analysis.
#' 
#' Input can be via environmental variables, or through direct
#' definition using `...` function arguments. In this case,
#' the direct function argument is used when available, and
#' environmental variable is used otherwise.
#' For example:
#' `get_slicejam_args(FC_FILE="test.fc")` will define `FC_FILE="test.fc"`
#' and will not use `Sys.getenv("FC_FILE")`.
#' 
#' Recognized arguments:
#' 
#' * FC_FILE - path to featureCounts .fc file
#' * CURATION_TXT - path to tab-delimited curation.txt file
#' * GTF - path to GTF file
#' * GTFNAME - optional short GTF name
#' * DRYRUN - 1/0 whether to enable dry-run mode
#' * OUTDIR - optional output directory
#' * ATAC - 1/0 whether to enable ATAC mode
#' * MGM - numeric max group mean threshold
#' * NORM - normalization method: quantile, median, mediangroup, none
#' * NORMMIN - normalization minimum_mean used for median, mediangroup
#' * GROUPCHECK - 1/0 whether to stop analysis after groups and contrasts
#' * SAVE_RDATA - 1/0 whether to save RData during processing
#' * UPSTREAM_PROMOTER - default 2000
#' * DOWNSTREAM_PROMOTER - default 200
#' * UPSTREAM_TTS - default 1000
#' * DOWNSTREAM_TTS - default 1000
#' * DETECTED_TX - optional path to detected transcripts, one per line
#' * DETECTED_GENES - optional path to detected genes, one per line
#' 
#' Returned arguments derived from input:
#' 
#' * fc_file - valid fc_file
#' * fc_filepath - valid full path to featureCounts file
#' * NORMSHORT - short normalization name
#' * fc_base - base filename derived from featureCounts filename
#' * fc_basedir - base directory derived from featureCounts filename
#' * fc_html - output HTML filename
#' * knit_root_dir - full path to directory used for knitr output
#' * cache_dir - directory used for knitr file cache
#' 
#' @export
get_slicejam_args <- function
(use_first_fc=TRUE,
 verbose=FALSE,
 type=c("printDebug",
    "data.frame"),
 ...)
{
   type <- match.arg(type);
   # FC_FILE - featureCounts file usually with .fc file extension
   inplist <- list(...);
   if ("FC_FILE" %in% names(inplist)) {
      fc_env <- inplist$FC_FILE
   } else {
      fc_env <- Sys.getenv("FC_FILE");
   }
   if (nchar(fc_env) > 0) {
      if (!file.exists(fc_env)) {
         stop(paste0("FC_FILE does not exist as a file:", fc_env,
            " in working directory:", getwd()));
      }
      fc_file <- fc_env;
   } else {
      if (use_first_fc) {
         fc_file <- head(list.files(pattern="[.]fc$"), 1);
      }
      if (length(fc_file) == 0) {
         stop(paste0(
            "No featureCounts .fc file was found in the working directory:",
            getwd()));
      }
   }
   
   # CURATION_TXT - curation tab-delimited text file
   # - to convert fc colnames to sample, group, batch
   if ("CURATION_TXT" %in% names(inplist)) {
      CURATION_TXT <- inplist$CURATION_TXT;
   } else {
      CURATION_TXT <- Sys.getenv("CURATION_TXT",
         "curation.txt");
   }
   if (length(CURATION_TXT) == 0 ||
         nchar(CURATION_TXT) == 0 ||
         !file.exists(CURATION_TXT)) {
      jamba::printDebug("Error: ",
         c("No ",
            "CURATION_TXT",
            " file was found in working directory:"),
         getwd(),
         c(" using ", "CURATION_TXT:"),
         CURATION_TXT,
         sep="",
         fgText=c("red", "dodgerblue"));
      CURATION_TXT <- normalizePath(CURATION_TXT);
      stop(paste0("No CURATION_TXT file was found in working directory:",
         getwd(), ", using CURATION_TXT:", CURATION_TXT));
   }
   
   # GTF used to annotate peaks for neighboring genes
   # GTFNAME optional short name for GTF file
   if ("GTF" %in% names(inplist)) {
      GTF <- inplist$GTF;
   } else {
      GTF <- Sys.getenv("GTF");
   }
   if ("GTFNAME" %in% names(inplist)) {
      GTFNAME <- inplist$GTFNAME;
   } else {
      GTFNAME <- Sys.getenv("GTFNAME");
   }
   if (length(GTF) > 0 && nchar(GTF) > 0) {
      if (!file.exists(GTF)) {
         stop(paste0("GTF file not found:", GTF));
      }
      if (length(GTFNAME) == 0 || nchar(GTFNAME) == 0) {
         GTFNAME <- paste0("_",
            gsub("gene[s]*", "",
               gsub("[-_ ]+", "",
                  gsub("[.][^.]+$", "", basename(GTF)))));
      }
      gtf_stem <- paste0("_gtf", GTFNAME);
      Sys.setenv(GTFNAME=GTFNAME);
   } else {
      stop("GTF was not supplied and is required.");
   }
   
   # DRYRUN not likely to be used here but defines "dry-run"
   # - which does not process data, only prints debug output
   if ("DRYRUN" %in% names(inplist)) {
      DRYRUN <- inplist$DRYRUN;
   } else {
      DRYRUN <- Sys.getenv("DRYRUN", "1");
   }
   if (!"0" %in% DRYRUN) {
      DRYRUN <- TRUE;
   } else {
      DRYRUN <- FALSE;
   }
   
   # GROUPCHECK whether to stop analysis after groups and contrasts
   if ("GROUPCHECK" %in% names(inplist)) {
      GROUPCHECK <- inplist$GROUPCHECK;
   } else {
      GROUPCHECK <- Sys.getenv("GROUPCHECK");
   }
   if (nchar(GROUPCHECK) > 0 && any(c("1", "true", "yes", "t", "y") %in% tolower(GROUPCHECK))) {
      GROUPCHECK <- TRUE;
   } else {
      GROUPCHECK <- FALSE;
   }

   # SAVE_RDATA
   if ("SAVE_RDATA" %in% names(inplist)) {
      SAVE_RDATA <- inplist$SAVE_RDATA;
   } else {
      SAVE_RDATA <- Sys.getenv("SAVE_RDATA");
   }
   if (length(SAVE_RDATA) == 0 || nchar(SAVE_RDATA) == 0 || any(c("0", "false", "f", "no", "n") %in% tolower(SAVE_RDATA))) {
      SAVE_RDATA <- FALSE;
   } else {
      SAVE_RDATA <- TRUE;
   }
   
   # DETECTED_TX, DETECTED_GENES
   if ("DETECTED_TX" %in% names(inplist)) {
      detectedTx_file <- inplist$DETECTED_TX;
   } else {
      detectedTx_file <- Sys.getenv("DETECTED_TX");
   }
   if (length(detectedTx_file) > 0 && nchar(detectedTx_file) > 0 && file.exists(detectedTx_file)) {
      detectedTx <- readLines(detectedTx_file);
   } else {
      detectedTx <- NULL;
   }
   if ("DETECTED_TX" %in% names(inplist)) {
      detectedGenes_file <- inplist$DETECTED_GENES;
   } else {
      detectedGenes_file <- Sys.getenv("DETECTED_GENES");
   }
   if (length(detectedGenes_file) > 0 && nchar(detectedGenes_file) > 0 && file.exists(detectedGenes_file)) {
      detectedGenes <- readLines(detectedGenes_file);
   } else {
      detectedGenes <- NULL;
   }
   

   # upstream and downstream ranges
   if ("UPSTREAM_PROMOTER" %in% names(inplist)) {
      UPSTREAM_PROMOTER <- inplist$UPSTREAM_PROMOTER;
   } else {
      UPSTREAM_PROMOTER <- Sys.getenv("UPSTREAM_PROMOTER", 2000);
   }
   if ("DOWNSTREAM_PROMOTER" %in% names(inplist)) {
      DOWNSTREAM_PROMOTER <- inplist$DOWNSTREAM_PROMOTER;
   } else {
      DOWNSTREAM_PROMOTER <- Sys.getenv("DOWNSTREAM_PROMOTER", 200);
   }
   if ("UPSTREAM_TTS" %in% names(inplist)) {
      UPSTREAM_TTS <- inplist$UPSTREAM_TTS;
   } else {
      UPSTREAM_TTS <- Sys.getenv("UPSTREAM_TTS", 1000);
   }
   if ("DOWNSTREAM_TTS" %in% names(inplist)) {
      DOWNSTREAM_TTS <- inplist$DOWNSTREAM_TTS;
   } else {
      DOWNSTREAM_TTS <- Sys.getenv("DOWNSTREAM_TTS", 1000);
   }
   upstream_promoter <- as.numeric(UPSTREAM_PROMOTER);
   downstream_promoter <- as.numeric(DOWNSTREAM_PROMOTER);
   upstream_tts <- as.numeric(UPSTREAM_TTS);
   downstream_tts <- as.numeric(DOWNSTREAM_TTS);
   
   # OUTDIR directory for Rmarkdown output
   if ("OUTDIR" %in% names(inplist)) {
      OUTDIR <- inplist$OUTDIR;
   } else {
      OUTDIR <- Sys.getenv("OUTDIR");
   }
   
   # CACHEDIR
   if ("CACHEDIR" %in% names(inplist)) {
      CACHEDIR <- inplist$CACHEDIR;
   } else {
      CACHEDIR <- Sys.getenv("CACHEDIR");
   }
   if (length(CACHEDIR) == 0 || nchar(CACHEDIR) == 0) {
      CACHEDIR <- file.path(OUTDIR, "cache", "");
   }
   if (!grepl("/$", CACHEDIR)) {
      CACHEDIR <- file.path(CACHEDIR, "");
   }
   
   # ATAC enables ATAC mode
   # - which requires genes to have an ATAC peak near the promoter
   #   during gene annotation of peaks
   if ("ATAC" %in% names(inplist)) {
      ATAC <- inplist$ATAC;
   } else {
      ATAC <- Sys.getenv("ATAC");
   }
   if (length(ATAC) == 0 || nchar(ATAC) == 0 || ATAC %in% c("0")) {
      ATAC <- 0;
   } else {
      ATAC <- 1;
   }
   
   # MGM maxGroupMean
   if ("MGM" %in% names(inplist)) {
      MGM <- inplist$MGM;
   } else {
      MGM <- Sys.getenv("MGM", 4);
   }
   MGM <- jamba::rmNA(as.numeric(unlist(strsplit(MGM, "[ ,]+"))));
   if (length(MGM) == 0) {
      MGM <- c(4);
   }
   
   # NORM normalization method
   if ("NORM" %in% names(inplist)) {
      NORM <- inplist$NORM;
   } else {
      NORM <- Sys.getenv("NORM");
   }
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

   # NORMMIN minimum_mean expression
   # - used when NORM includes median normalization
   if ("NORMMIN" %in% names(inplist)) {
      NORMMIN <- inplist$NORMMIN;
   } else {
      NORMMIN <- Sys.getenv("NORMMIN");
   }
   if (nchar(NORMMIN) == 0) {
      NORMMIN <- MGM;
   }
   
   # Define base filename used for output
   fc_base <- gsub("[.]fc$|[.]fc[.]txt$", "",
      fc_file);
   fc_filepath <- normalizePath(fc_file);
   fc_dirname <- dirname(fc_filepath);
   
   ## add the max group mean threshold (mgm) _mgm4
   fc_base <- paste0(fc_base, "_mgm", MGM);
   ## add the normalization short name
   fc_base <- paste0(fc_base, "_", NORMSHORT);
   ## add the ATAC mode, _atac1 when ATAC-mode is enabled
   if (!ATAC %in% c(0)) {
      fc_base <- paste0(fc_base, "_atac", ATAC);
   }
   # add the GTF stem to the base name
   fc_base <- paste0(fc_base, gtf_stem);
   
   # knitr root path auto-generated or using OUTDIR directly
   # - either way it uses a fully described file path
   if (nchar(OUTDIR) == 0) {
      fc_basedir <- paste0(fc_base,  "_analysis");
      knit_root_dir <- file.path(fc_dirname,
         fc_basedir);
      # knitr root directory
      OUTDIR <- knit_root_dir;
   } else {
      OUTDIR <- normalizePath(OUTDIR,
         mustWork=FALSE);
      fc_base <- OUTDIR;
      fc_basedir <- OUTDIR;
      # knitr root directory
      knit_root_dir <- OUTDIR;
   }
   
   # html output filename
   fc_html <- file.path(knit_root_dir,
      paste0(basename(fc_base),  ".html"));
   # cache directory
   CACHEDIR <- file.path(knit_root_dir, "cache", "");
   
   sliceargs <- list();
   sliceargs$fc_file <- fc_file;
   sliceargs$fc_filepath <- fc_filepath;
   sliceargs$CURATION_TXT <- CURATION_TXT;
   sliceargs$GTF <- GTF;
   sliceargs$GTFNAME <- GTFNAME;
   sliceargs$DRYRUN <- DRYRUN;
   sliceargs$OUTDIR <- OUTDIR;
   sliceargs$CACHEDIR <- CACHEDIR;
   sliceargs$ATAC <- ATAC;
   sliceargs$MGM <- MGM;
   sliceargs$NORM <- NORM;
   sliceargs$NORMMIN <- NORMMIN;
   sliceargs$NORMSHORT <- NORMSHORT;
   sliceargs$GROUPCHECK <- GROUPCHECK;
   sliceargs$SAVE_RDATA <- SAVE_RDATA;
   
   sliceargs$UPSTREAM_PROMOTER <- upstream_promoter;
   sliceargs$DOWNSTREAM_PROMOTER <- downstream_promoter;
   sliceargs$UPSTREAM_TTS <- upstream_tts;
   sliceargs$DOWNSTREAM_TTS <- downstream_tts;
   
   sliceargs$DETECTED_TX <- detectedTx;
   sliceargs$DETECTED_GENES <- detectedGenes;
   
   sliceargs$fc_base <- fc_base;
   sliceargs$fc_basedir <- fc_basedir;
   sliceargs$fc_html <- fc_html;
   sliceargs$knit_root_dir <- knit_root_dir;
   if (verbose) {
      debug_df <- printDebugList(sliceargs,
         type=type);
      if ("data.frame" %in% type) {
         # print kable
         print(knitr::kable(
            debug_df,
            caption=paste0(
               "slicejam arguments:")));
         
      }
   }
   return(sliceargs);
}

#' Print debug with list input
#' 
#' Print debug with list input
#' 
#' This function is a lightweight wrapper around
#' `jamba::printDebug()` that allows list input,
#' in a way the prints each list name.
#' 
#' @param ... items which will be printed
#' @param maxprint `integer` maximum number of elements to print
#'    for each item in `...`.
#' @param type `character` string indicating the type of output:
#'    * `"printDebug"` - calls `jamba::printDebug()` to print to the R console
#'    * `"data.frame"` - assembles a `data.frame`
#' 
#' @examples
#' printDebugList(A=list(B=1:3, C=4:9, D=letters[1:3]),
#'    something_else=sample(LETTERS, 5))
#'
#' printDebugList(A=list(B=1:3, C=4:9, D=letters[1:3]),
#'    something_else=sample(LETTERS, 5),
#'    maxprint=4)
#' 
#' printDebugList(A=list(B=1:3, C=4:9, D=letters[1:3]),
#'    E=list(1:5, G=1:4),
#'    something_else=sample(LETTERS, 5),
#'    maxprint=4)
#' 
#' @export
printDebugList <- function
(...,
 maxprint=10,
 type=c("printDebug",
    "data.frame"),
 debug=FALSE)
{
   type <- match.arg(type);
   
   # determine which arguments are intended for jamba::printDebug()
   # all others will be printed onscreen
   fnames <- setdiff(names(formals(jamba::printDebug)), "...");
   inplist <- list(...);
   if (length(names(inplist)) == 0) {
      printlist <- inplist;
      arglist <- NULL;
   } else {
      printlist <- inplist[!names(inplist) %in% fnames];
      arglist <- inplist[names(inplist) %in% fnames];
   }
   
   indent <- "";
   
   ############################################
   # custom print function for nested lists
   ret_df <- data.frame(argument=NA_character_,
      value=NA_character_)[0,,drop=FALSE];
   print_list_arg <- function
   (printlist,
    i,
    indent="",
    type=c("printDebug",
       "data.frame"))
   {
      type <- match.arg(type);
      if (length(names(printlist)) > 0 && nchar(names(printlist)[i]) > 0) {
         pname <- paste0(indent,
            names(printlist)[i],
            ":");
      } else {
         pname <- "[list]:";
      }
      ret_df <- data.frame(argument=NA_character_,
         value=NA_character_)[0,,drop=FALSE];
      printvals <- NULL;
      if (debug) {
         print("1   arglist:");print(arglist);
         print("1     pname:");print(pname);
         print("1 printvals:");print(printvals);
      }
      if ("printDebug" %in% type) {
         do.call(jamba::printDebug,
            c(arglist,
               list(pname)));
      } else if ("data.frame" %in% type) {
         ret_df <- rbind(ret_df,
            data.frame(
               argument=jamba::rmNULL(pname,
                  nullValue=""),
               value=jamba::cPaste(printvals, sep=", ")));
      }
      
      indent <- paste0("   ", indent);
      
      for (j in seq_along(printlist[[i]])) {
         if (is.list(printlist[[i]][[j]])) {
            new_df <- print_list_arg(printlist[[i]], j, indent, type);
            if ("data.frame" %in% type) {
               ret_df <- rbind(ret_df,
                  new_df);
            }
         } else {
            if (length(names(printlist[[i]])) > 0 && nchar(names(printlist[[i]])[j]) > 0) {
               pname <- paste0(indent,
                  names(printlist[[i]])[j],
                  ":");
            } else {
               pname <- indent;
            }
            printvals <- printlist[[i]][[j]];
            if (length(printvals) > maxprint) {
               printvals <- c(head(printvals, maxprint),
                  paste0("...(", jamba::formatInt(length(printvals)), " total)"));
            }
            if (debug) {
               print("2   arglist:");print(arglist)
               print("2     pname:");print(pname)
               print("2 printvals:");print(printvals)
            }
            if ("printDebug" %in% type) {
               do.call(jamba::printDebug,
                  c(arglist,
                     list(pname),
                     list(printvals)));
            } else if ("data.frame" %in% type) {
               ret_df <- rbind(ret_df,
                  data.frame(
                     argument=jamba::rmNULL(pname,
                        nullValue=""),
                     value=jamba::cPaste(printvals, sep=", ")));
            }
         }
      }
      indent <- sub("   ", "", indent);
      return(ret_df);
   }
   ############################################
   
   for (i in seq_along(printlist)) {
      if (is.list(printlist[[i]])) {
         new_df <- print_list_arg(printlist, i, indent, type);
         if ("data.frame" %in% type) {
            ret_df <- rbind(ret_df,
               new_df);
         }
      } else {
         if (length(names(printlist)) > 0 && nchar(names(printlist)[i]) > 0) {
            pname <- paste0(names(printlist)[i], ":");
         } else {
            pname <- NULL;
         }
         printvals <- printlist[[i]];
         #if (length(printvals) > maxprint) {
         #   printvals <- c(head(printvals, maxprint),
         #      paste0("... (", jamba::formatInt(length(printvals)), " total)"));
         #}
         if (debug) {
            print("3   arglist:");print(arglist)
            print("3     pname:");print(pname)
            print("3 printvals:");print(printvals)
         }
         if ("printDebug" %in% type) {
            do.call(jamba::printDebug,
               c(arglist,
                  list(pname), list(printvals)));
         } else if ("data.frame" %in% type) {
            ret_df <- rbind(ret_df,
               data.frame(
                  argument=jamba::rmNULL(pname,
                     nullValue=""),
                  value=jamba::cPaste(printvals, sep=", ")));
         }
      }
   }
   if ("data.frame" %in% type) {
      return(ret_df)
   }
}

