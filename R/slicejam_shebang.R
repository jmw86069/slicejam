
#' Update slicejam Rscript shebang
#' 
#' Update slicejam Rscript shebang
#' 
#' This function updates the shebang header line in the Rscript
#' files installed with the slicejam R package.
#' 
#' The shebang makes the Rscript file properly executable
#' on a command-line by pointing directly to the correct
#' `Rscript` utility.
#' 
#' Note this function is somewhat experimental, since the
#' process of updating an installed R package script with
#' a shebang line should be a formal part of the build
#' process, but this process is not well documented.
#' 
#' @family slicejam commandline
#' 
#' @export
update_slicejam_shebang <- function
(dryrun=FALSE,
 verbose=TRUE,
 ...)
{
   ## Direct R_HOME method to determine Rscript
   rscript <- paste0(
      Sys.getenv("R_HOME"),
      "/bin",
      Sys.getenv("R_ARCH_BIN"),
      "/Rscript");
   file.exists(rscript);
   rscript <- file.path(R.home(), "bin/Rscript");

   ## if it fails, fallback by other methods
   if (!file.exists(rscript)) {

      ## Find current R executable
      args_v <- commandArgs();
      jamba::printDebug("R command:\n",
         args_v[1]);
      
      if (grepl("RStudio", args_v[1], ignore.case=TRUE)) {
         ## When running inside Rstudio
         rscript <- head(
            list.files(pattern="^Rscript$",
               full.names=TRUE,
               path=c(R.home(), list.dirs(R.home(), recursive=FALSE))),
            1);
      } else {
         ## When running inside R (not RStudio)
         rpaths <- c(dirname(dirname(args_v[1])),
            dirname(args_v[1]));
         rscript <- list.files(path=rpaths,
            pattern="[rR]script",
            full.names=TRUE);
      }
   }
   
   if (length(rscript) == 0 || nchar(rscript) == 0 || !file.exists(rscript)) {
      stop("Did not find Rscript.");
   }
   
   ## new shebang
   new_shebang <- paste0("#!", rscript);
   if (verbose) {
      jamba::printDebug("update_slicejam_shebang(): ",
         "Rscript shebang:\n",
         new_shebang);
   }

   ## update these files with the specific path to Rscript
   scripts_to_update <- list.files(
      path=system.file(package="slicejam", "exec", "bin"),
      full.names=TRUE);
   if (verbose) {
      jamba::printDebug("update_slicejam_shebang(): ",
         "scripts_to_update:\n",
         scripts_to_update,
         sep="\n");
   }

   for (iscript in scripts_to_update) {
      script_lines <- readLines(iscript);
      if (verbose) {
         jamba::printDebug("update_slicejam_shebang(): ",
            "iscript:\n",
            iscript);
      }
      if (grepl("^[#][#!]", script_lines[1])) {
         if (verbose) {
            jamba::printDebug("update_slicejam_shebang(): ",
               "   the shebang line was updated from:\n",
               script_lines[1]);
         }
         script_lines[1] <- new_shebang;
         if (verbose) {
            jamba::printDebug("update_slicejam_shebang(): ",
               "   the shebang line was updated to:\n",
               script_lines[1]);
         }
         if (dryrun) {
            cat(paste(script_lines[1:4], collapse="\n"), "\n");
         }
         #writeLines(con=iscript,
         #   text=script_lines);
      }
   }
   
}
