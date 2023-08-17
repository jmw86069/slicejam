
#' Setup slicejam analysis scripts
#' 
#' Setup slicejam scripts for command-line analysis
#' 
#' @family slicejam commandline
#' 
#' @param output_path `character` file path to a directory where the
#'    resulting `.R` and `.Rmd` files should be stored.
#' @param overwrite `logical` indicating whether to overwrite existing
#'    files if they already exist in the `output_path` directory.
#' @param dryrun `logical` indicating whether to print commands without
#'    performing the file copy steps.
#' @param bash_shebang `character` line used at the top of the output
#'    script, traditionally called a "shebang". This line points
#'    the script to a specific executable environment, by default
#'    points to bash. This argument is only appropriate for linux- and
#'    Unix-like operating systems.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @returns `character` string invisibly returned, which contains
#'    the full bash shell script.
#' 
#' @examples
#' if (FALSE) {
#'    testdir <- file.path(tempdir(), "slicejam_testing")
#'    if (!dir.exists(testdir)) {
#'       dir.create(testdir)
#'    }
#'    # print steps to be taken without performing them
#'    setup_slicejam(testdir, dryrun=TRUE)
#'    
#'    # performing the steps
#'    setup_slicejam(testdir, dryrun=FALSE)
#' }
#' 
#' @export
setup_slicejam <- function
(output_path=".",
 overwrite=FALSE,
 dryrun=TRUE,
 bash_shebang="#!/usr/bin/env bash",
 wrapper_script="run_slicejam.sh",
 verbose=TRUE,
 ...)
{
   # validate output_path
   if (length(output_path) == 0 || nchar(head(output_path, 1)) == 0) {
      output_path <- ".";
   } else {
      output_path <- head(output_path, 1)
   }
   normalized_output_path <- normalizePath(output_path);
   if (verbose) {
      jamba::printDebug("", #"setup_slicejam(): ",
         "output_path:            '",
         output_path, "'")
      jamba::printDebug("", #"setup_slicejam(): ",
         "normalized_output_path: '",
         normalized_output_path, "'")
   }

   if (length(wrapper_script) == 0 || nchar(wrapper_script[1]) == 0) {
      wrapper_script <- "run_slicejam.sh"
   }
   output_wrapper <- file.path(normalized_output_path,
      wrapper_script)
   
   # slicejam package dir
   slicejam_exec <- system.file(package="slicejam", "exec");
   if (verbose) {
      jamba::printDebug("", #"setup_slicejam(): ",
         "slicejam exec path:     '",
         slicejam_exec, "'")
   }
   
   # determine RHOME
   RHOME <- Sys.getenv("R_HOME")
   
   # define R_LIBS
   # TODO: delimiter should be semicolon on Windows
   R_LIBS <- jamba::cPaste(.libPaths(),
      sep=":");

   ######################################
   # helper function: check output file
   #
   # returns TRUE if file can be copied
   check_output_file <- function
   (output_file,
    overwrite=FALSE,
    dryrun=FALSE,
    verbose=FALSE)
   {
      #
      copy_filename <- basename(output_file)
      if (file.exists(output_file)) {
         if (TRUE %in% overwrite) {
            if (verbose) {
               jamba::printDebug("", #"setup_slicejam(): ",
                  "overwriting: ",
                  copy_filename)
            }
         } else {
            if (verbose) {
               jamba::printDebug("", #"setup_slicejam(): ",
                  "skipping: ",
                  copy_filename,
                  fgText=c("darkorange1",
                     "dodgerblue",
                     "firebrick3"))
            }
            return(FALSE);
         }
      } else {
         if (verbose) {
            jamba::printDebug("", #"setup_slicejam(): ",
               "copying: ",
               copy_filename)
         }
      }
      if (TRUE %in% dryrun) {
         if (verbose) {
            jamba::printDebug("", #"setup_slicejam(): ",
               "skipped due to dryrun=TRUE.",
               indent=6,
               fgText=c("darkorange2", "firebrick3"))
         }
         return(FALSE)
      }
      return(TRUE)
   }
   
   # copy run_slicejam.R
   copy_filenames <- c("run_slicejam.R",
      "slicejam_analysis.Rmd")
   for (copy_filename in copy_filenames) {
      output_file <- file.path(output_path,
         copy_filename)
      source_file <- file.path(slicejam_exec,
         copy_filename)
      # check for output file handling
      do_copy <- check_output_file(output_file,
         overwrite=overwrite,
         dryrun=dryrun,
         verbose=verbose)
      if (TRUE %in% do_copy) {
         file.copy(from=source_file,
            to=output_file,
            overwrite=TRUE)
      }
   }
   
   # define absolute path to run_slicejam.R
   run_slicejam_path <- file.path(
      normalized_output_path,
      "run_slicejam.R")
   
   # now create commandline script
   bash_script <- paste0(
bash_shebang, '
# the line above should contain shebang which points to bash
b="\033[1m";
r="\033[0m";
h="\033[1m\033[36m";
u="\033[1m\033[37m";
export R_LIBS="', R_LIBS, '";
export RHOME="', RHOME, '";
DRYRUN=${DRYRUN:-"1"};
GROUPCHECK=${GROUPCHECK:-"0"};

function usage {
   echo "${h}run_slicejam.sh ${u}version ', packageVersion("slicejam"), '${r}";
   echo "   Wrapper utility to perform slicejam analysis.";
   echo "   Settings are defined as environment variables.";
   echo "";
   echo "${h}Processing options${r}:";
   echo "   ${b}DRYRUN${r}: ${u}0/1${r} indicating whether to perform dry-run analysis";
   echo "   ${b}GROUPCHECK${r}: ${u}0/1${r} whether to stop at sample group information";
   echo "${h}Data input options${r}:";
   echo "   ${b}FC_FILE${r}: ${u}path${r} to tab-delimited text featureCounts file,";
   echo "      expected to have file extension ${u}.fc${r} or ${u}.fc.txt${r}";
   echo "   ${b}GTF${r}: ${u}path${r} to gene GTF file";
   echo "   ${b}GTFNAME${r}: optional label to use for the GTF file";
   echo "   ${b}CURATION_TXT${r}: ${u}path${r} to tab-delimited sample curation,";
   echo "      with column names: ${u}Filename, Sample, Run, Group${r}."
   echo "      The Filename values should match a substring of the";
   echo "      column headers in the featureCounts file.";
   echo "   ${b}MASK${r}: ${u}path${r} (optional) to BED file with regions to flag as masked.";
   echo "   ${b}DETECTED_TX${r}: ${u}path${r} (optional) to detected transcripts,";
   echo "      one per line, matched with transcript_id entries in the GTF.";
   echo "   ${b}DETECTED_GENES${r}: ${u}path${r} (optional) to detected genes,";
   echo "      one per line, matched with gene_name entries in the GTF.";
   echo "${h}Analysis options${r}:";
   echo "   ${b}ATAC${r}: ${u}0/1${r} whether to enable ATAC-mode promoter annotation";
   echo "   ${b}MGM${r}: ${u}numeric${r} max group mean noise threshold";
   echo "      during processing, default ${u}MGM=4${r}";
   echo "   ${b}NORM${r}: ${u}name${r} of normalization method to use${r}:";
   echo "      ${u}quantile${r} (default)";
   echo "      ${u}median${r}: median log-ratio normalization";
   echo "      ${u}mediangroup${r}: median log-ratio normalization within group only";
   echo "      ${u}none${r}: no additional normalization";
   echo "   ${b}NORMMIN${r}: ${u}numeric${r} (optional) mean signal threshold";
   echo "      during normalization by median and mediangroup, default is ${u}MGM${r}";
   echo "   ${b}UPSTREAM_PROMOTER,DOWNSTREAM_PROMOTER${r}: ${u}numeric${r} ";
   echo "      promoter range relative to the TSS, default ${u}2000,200${r}";
   echo "   ${b}UPSTREAM_TTS,DOWNSTREAM_TTS${r}: ${u}numeric${r}";
   echo "      range relative to the TTS, default ${u}1000,1000${r}";
   echo "${h}Output options${r}:";
   echo "   ${b}OUTDIR${r}: ${u}path${r} (optional) custom name for the";
   echo "   output directory.";
   echo "   ${b}SAVE_RDATA${r}: ${u}0/1${r} whether to save .RData file for re-use.";
   echo "      default ${u}0${r} since cache is already saved.";
   echo "";
   echo "${h}Example dry-run${r}:
   ${u}DRYRUN=1 GROUPCHECK=1 FC_FILE=some_counts_file.fc GTF=gencode_v32.gtf \\\\
   CURATION_TXT=curation.txt OUTDIR=slicejam_output ./run_slicejam.sh${r}";
   echo "";
   echo "${h}Analysis that stops at sample grouping${r}:
   ${u}DRYRUN=0 GROUPCHECK=1 FC_FILE=some_counts_file.fc GTF=gencode_v32.gtf \\\\
   CURATION_TXT=curation.txt OUTDIR=slicejam_output ./run_slicejam.sh${r}";
   echo "";
   echo "${h}Full analysis${r}:
   ${u}DRYRUN=0 GROUPCHECK=0 FC_FILE=some_counts_file.fc GTF=gencode_v32.gtf \\\\
   CURATION_TXT=curation.txt OUTDIR=slicejam_output ./run_slicejam.sh${r}";
   exit;
}

if [[ ( "--help" == "${1}" || "-h" == "${1}" ) ]]; then
   usage;
fi;
if [[ ( "" == "${FC_FILE}" || ! -f "${FC_FILE}" ) ]]; then
   echo "Error: No FC_FILE was defined.";
   usage;
fi;

DRYRUN="${DRYRUN}" \\
   GROUPCHECK="${GROUPCHECK}" \\
   GTF="${GTF}" \\
   FC_FILE="${FC_FILE}" \\
   CURATION_TXT="${CURATION_TXT}" \\
   OUTDIR="${OUTDIR}" \\
   Rscript ',
      run_slicejam_path, '
');

   # check for output file handling
   do_copy <- check_output_file(output_wrapper,
      overwrite=overwrite,
      dryrun=dryrun,
      verbose=verbose)
   if (TRUE %in% do_copy) {
      cat(file=output_wrapper,
         bash_script,
         sep="\n");
   }

   return(invisible(bash_script));
}
