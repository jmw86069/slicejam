
#' Convert featureCount .fc file to curation_txt data.frame
#' 
#' Convert featureCount .fc file to curation_txt data.frame, used to
#' populate reasonably close approximation of a sample curation.txt file.
#' 
#' This function is intended to be a helper function to create an
#' empty `"curation.txt"` file to use in slicejam. The output should
#' contain these columns at minimum:
#' 
#' * `"Filename"` - the first column should contain the base filename,
#' without the directory (path) to the file, only the filename.
#' * `"Run"` - optional, but used when multiple sequencing runs, and/or
#' sequencing machines were used to produce the sequencing data.
#' * `"Group"` - column containing `character` group names, typically
#' in the form where each experiment factor is separated by underscore `"_"`,
#' For example `"Dex_KO"` and `"Veh_WT"`. these values are used in slicejam
#' to intuit appropriate statistical contrasts by comparing changes
#' in only one factor at a time.
#' * `"Label"` - column containing unique `character` values for each row,
#' used to define a label for each sample, and is displayed on each
#' figure.
#' 
#' The procedure used to create the `data.frame`:
#' 
#' * The directory is removed from each filename, by calling `basename()`.
#' * The filename extension is removed by matching `suffix_extension`
#' and removing all text at, and to the right side of the match.
#' * Filenames are edited by substituting `curate_from` with `curate_to`,
#' if defined.
#' * Filenames have any non-alphanumeric or underscore "_" to `"."`.
#' * Remaining Filenames are split by argument `split` whose default
#' is `"[-_. ]+"`. The resulting fields are assembled into a `data.frame`.
#' * Columns are scanned using `run_prefix` to determine which column
#' (if any) has the most matches to `run_prefix`. The column with the
#' highest number of matching fields is renamed `"Run"`.
#' * The same procedure is used with `rep_prefix`, and the winning column
#' is renamed `"Run"`.
#' * Columns with only one value are removed from group colnames, and remaining
#' columns are concatenated with delimiter `"_"` into a group name.
#' These columns are removed, and replaced with one column `"Group"` with
#' these group values. If no columns survived to this point, a column
#' `"Group"` is created with empty string `""`, to be edited by the user.
#' * If `"Group"` is non-empty, it is used to create a `"Label"` column:
#' 
#'    1. If `"Rep"` exists, it is appended to the group value.
#'    2. If `"Run"` exists, it is also appended to the group.
#'    3. The combined values are forced unique by `jamba::makeNames()`
#'    which by default uses suffix `"_v#"` to any non-unique labels.
#' 
#' @returns `data.frame` with colnames suitable for use as argument
#'    `curation_txt` in `"slicejam_analysis.Rmd"`. The `data.frame` should
#'    be edited as relevant for downstream analysis.
#' 
#' @family slicejam utilities
#'
#' @param fcfile `character` with either:
#'    1. path to a featureCounts tab-delimited text file, OR
#'    2. vector with column headers directly.
#' @param fc_colnames `character` vector of featureCounts column headers
#'    used to indicate genome coordinates, which are not sample headers.
#' @param suffix_extension `character` vector of file extension patterns,
#'    used to trim the right side of each filename, removing anything that
#'    matches ".extension" and everything to the right of this pattern.
#'    The "." is added automatically.
#' @param curate_from vector of regular expression
#'    patterns used to edit the trimmed filename, intended to correct
#'    mistakes in the filename, or to edit potential delimiters before
#'    the fields are split into multiple columns of a `data.frame`.
#' @param curate_to `character` vector of values used for replacement for
#'    each matching value in `curate_from`.
#' @param run_prefix `character` prefix used to match the beginning of each
#'    string in each column of the `data.frame`, used to recognize
#'    the most likely column that contains sequencing run values.
#'    The prefix values should be something like the machine project number,
#'    in the form `"NOV0399"`, `"NOVA0399"`, `"NS52099"`, `"MS0032"`.
#' @param rep_prefix `character` prefix used to match the beginning of each
#'    string in each column of the `data.frame`, used to recognize
#'    the most likely column that contains replicate number values.
#'    The prefix values should be something like "rep", or "R",
#'    to match values in the form `"rep1"`, `"Rep2"`, `"R1"`, `"v2"`.
#' @param split `character` string used by `strsplit()` to split filenames
#'    into multiple fields, which are then converted to a `data.frame`.
#' @param curation_out_file `character` path to output file, or `NULL`
#'    by default not to save a file. When defined, the output is saved
#'    to this file as tab-delimited text.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' fcfile <- path.expand("../112822_KO_ATAC_cov/UL3_KO_ATAC_1hDex_vs_EtOH_1000s_mergeBed_36files_gap25_slop0_max1000_num1_28nov2022.fc");
#' 
#' fcfile <- c("/some/file/path/NOVA0021_dH1_cloneA_DEX_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneA_DEX_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneA_VEH_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneA_VEH_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneB_DEX_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneB_DEX_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneB_DEX_R3_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneB_DEX_R4_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneB_VEH_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneB_VEH_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneB_VEH_R3_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_dH1_cloneB_VEH_R4_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_UL3_DEX_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_UL3_DEX_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_UL3_VEH_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NOVA0021_UL3_VEH_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50462_dH1_cloneA_VEH_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50462_dH1_cloneA_VEH_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50462_dH1_cloneB_VEH_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50462_dH1_cloneB_VEH_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50462_UL3_VEH_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50462_UL3_VEH_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneA_DEX_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneA_DEX_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneA_VEH_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneA_VEH_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneA_VEH_R3._ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneB_DEX_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneB_DEX_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneB_DEX_R3._ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneB_VEH_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_dH1_cloneB_VEH_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_UL3_DEX_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_UL3_DEX_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_UL3_VEH_R1_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam",
#'    "/some/file/path/NS50644_UL3_VEH_R2_ATAC.trimmed.clean.optdedupe.pcrdedupe.T2T.mapped.sorted.bam")
#' fc_to_curation(fcfile)
#' 
#' curate_from <- c("dH1_clone");
#' curate_to <- c("dH1");
#' fc_to_curation(fcfile, curate_from=curate_from, curate_to=curate_to)
#' 
#' @export
fc_to_curation <- function
(fcfile,
 fc_colnames=c("Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length"),
 suffix_extensions=c("trim",
    "mapped",
    "aligned",
    "sort",
    "bam",
    "sam",
    "fastq",
    "clean",
    "dedupe"),
 curate_from=c(),
 curate_to=c(),
 run_prefix=c("NOV",
    "NOVA",
    "NS",
    "MS",
    "MiSeq",
    "HS",
    "HiSeq",
    "NextSeq",
    "NovaSeq"),
 rep_prefix=c("R",
    "Rep",
    "v"),
 split="[-_. ]+",
 curation_out_file=NULL,
 ...)
{
   #
   if (length(fcfile) == 1) {
      df <- data.table::fread(fcfile, nrows=10, data.table=FALSE)
      data_colnames <- setdiff(colnames(df), fc_colnames)
   } else {
      data_colnames <- fcfile;
   }

   # trim suffix_extensions
   gsub_pattern <- paste0("[.](",
      jamba::cPaste(suffix_extensions, sep="|"),
      ").*$");
   sample_colnames <- gsub(gsub_pattern, "", basename(data_colnames))
   # Convert non-alphanumeric characters to "."
   filenames <- gsub("[^a-zA-Z0-9_]", ".", sample_colnames);
   
   # curate_from replacement with curate_to
   if (length(curate_from) > 0) {
      if (length(curate_to) == 0) {
         curate_to <- ""
      }
      curate_to <- rep(curate_to, length.out=length(curate_from))
      for (i in seq_along(curate_from)) {
         sample_colnames <- gsub(curate_from[i], curate_to[i], sample_colnames)
      }
   }

   # Convert non-alphanumeric characters to "."
   sample_colnames <- gsub("[^a-zA-Z0-9_]", ".", sample_colnames)
   
   # split into data.frame columns
   split_df <- data.frame(stringsAsFactors=FALSE,
      jamba::rbindList(strsplit(sample_colnames, "[-_. ]+")))
   split_colnames <- colnames(split_df);
   
   # detect best guess "Run" column
   run_grep <- paste0("^(", jamba::cPaste(run_prefix, sep="|"), ")[0-9]+");
   run_grep_count <- sapply(split_colnames, function(icol){
      sum(grepl(run_grep, ignore.case=TRUE, split_df[[icol]]))
   })
   run_colname <- NULL;
   if (any(run_grep_count > 0)) {
      run_colname <- names(which.max(run_grep_count));
      split_df <- jamba::renameColumn(split_df,
         from=run_colname,
         to="Run");
      split_colnames <- setdiff(split_colnames, run_colname);
   }
   
   # detect best guess "Run" column
   rep_grep <- paste0("^(", jamba::cPaste(rep_prefix, sep="|"), ")[0-9]+");
   rep_grep_count <- sapply(split_colnames, function(icol){
      sum(grepl(rep_grep, ignore.case=TRUE, split_df[[icol]]))
   })
   # define "Rep" if more than half the rows match the pattern
   rep_colname <- NULL;
   if (any(rep_grep_count > (nrow(split_df) / 2))) {
      rep_colname <- names(which.max(rep_grep_count));
      split_df <- jamba::renameColumn(split_df,
         from=rep_colname,
         to="Rep");
      split_colnames <- setdiff(split_colnames, rep_colname);
   }
   
   # trim (remove) columns with only one value
   if (length(split_colnames) > 0) {
      nvals <- sapply(split_colnames, function(icol){
         length(unique(split_df[[icol]]))
      })
      if (any(nvals == 1)) {
         drop_colnames <- names(nvals)[nvals <= 1]
         split_colnames <- setdiff(split_colnames,
            drop_colnames)
         # drop all but these columns
         keep_colnames <- setdiff(colnames(split_df),
            drop_colnames);
         # split_df <- split_df[, keep_colnames, drop=FALSE];
      }
   }
   
   # the rest are assumed to be group_colnames
   group_colname <- NULL;
   if (length(split_colnames) == 1) {
      split_df2 <- jamba::renameColumn(split_df,
         from=split_colnames,
         to="Group");
      group_colname <- "Group";
   } else if (length(split_colnames) > 1) {
      group_values <- jamba::pasteByRow(split_df[, split_colnames, drop=FALSE])
      split_df$Group <- group_values;
      keep_colnames <- setdiff(colnames(split_df),
         split_colnames);
      split_df <- split_df[, keep_colnames, drop=FALSE];
      group_colname <- "Group";
   } else {
      split_df$Group <- "";
   }
   
   # Create Label column
   if (length(group_colname) == 1) {
      label_colnames <- "Group"
      if (length(rep_colname) == 1) {
         label_colnames <- c(label_colnames, "Rep")
      }
      if (length(run_colname) == 1) {
         label_colnames <- c(label_colnames, "Run")
      }
      label_values <- jamba::makeNames(
         jamba::pasteByRow(split_df[, label_colnames, drop=FALSE]))
      split_df$Label <- label_values;
   }
   
   curation_txt <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      Filename=filenames,
      split_df)
   
   # optionally save output
   if (length(curation_out_file) == 1 && is.character(curation_out_file)) {
      tryCatch({
         data.table::fwrite(file=curation_out_file,
            x=curation_txt,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
      }, error=function(e){
         jamba::printDebug("fc_to_curation(): ",
            "File could not be saved.");
         print(e);
         NULL;
      })
   }
   
   curation_txt
}
