
## slicejam custom functions
## likely to be a placeholder until functions are moved to the
## appropriate jam R package

#' Curate vector to data.frame by pattern matching
#' 
#' Curate vector to data.frame by pattern matching
#' 
#' This function takes a vector, and a `data.frame` with one
#' column containing text string patterns. Whenever a pattern
#' matches the input vector, the annotations in `data.frame`
#' are used.
#' 
#' @examples
#' df <- data.frame(
#'    pattern=c("NOV14_p2w5_VEH",
#'       "NOV14_p4w4_VEH",
#'       "NOV14_UL3_VEH",
#'       "NS644_UL3VEH",
#'       "NS50644_UL3VEH",
#'       "NS644_p2w5VEH"),
#'    batch=c("NOV14",
#'       "NOV14",
#'       "NOV14",
#'       "NS644",
#'       "NS50644",
#'       "NS644"),
#'    group=c("p2w5_Veh",
#'       "p4w4_Veh",
#'       "UL3_Veh",
#'       "UL3_Veh",
#'       "UL3_Veh",
#'       "p2w5_Veh")
#' );
#' ## review the input table format
#' print(df);
#' x <- c("NOV14_p2w5_VEH_25_v2_CoordSort_deduplicated_SingleFrag_38to100.bam",
#'    "NOV14_p4w4_VEHrep1_25_v2_CoordSort_deduplicated_SingleFrag_38to100.bam",
#'    "NOV14_UL3_VEH_25_v2_CoordSort_deduplicated_SingleFrag_38to100.bam",
#'    "NS644_UL3VEH_25_v3_CoordSort_deduplicated_SingleFrag_38to100.bam",
#'    "NOV14_p2w5_VEH_50_v2_CoordSort_dedup_singleFragment.bam",
#'    "NOV14_UL3_VEH_50_v2_CoordSort_dedup_singleFragment.bam",
#'    "NS50644_UL3VEH_25_v3_CoordSort_deduplicated_SingleFrag.bam",
#'    "NS644_p2w5VEH_12p5_v3_CoordSort_deduplicated_SingleFrag_38to100.bam")
#' 
#' df_new <- curate_to_df_by_pattern(x, df);
#' ## Review the curated output
#' print(df_new);
#' 
#' ## Print a colorized image
#' colorSub <- colorjam::group2colors(unique(unlist(df_new)));
#' colorSub <- jamba::makeColorDarker(colorSub, darkFactor=-1.6, sFactor=-1.6);
#' k <- c(1,2,3,4,5,5,5,5);
#' df_colors <- as.matrix(df_new[,k]);
#' df_colors[] <- colorSub[df_colors];
#' opar <- par("mar"=c(3,3,4,3));
#' jamba::imageByColors(df_colors,
#'    adjustMargins=FALSE,
#'    cellnote=df_new[,k],
#'    flip="y",
#'    cexCellnote=c(0.4,0.5)[c(1,2,2,2,1,1,1,1)],
#'    xaxt="n",
#'    yaxt="n",
#'    groupBy="row");
#' axis(3,
#'    at=c(1,2,3,4,6.5),
#'    labels=colnames(df_new));
#' par(opar);
#' 
#' @param x `character` vector of input data, often filenames.
#' @param df `data.frame` containing a column with patterns, and
#'    one or more columns that contain annotations.
#' @param pattern_colname,group_colname,id_colname `character` string
#'    indicating colname to use for patterns, group, and identifier,
#'    respectively. The `group_colname` and `id_colname` may be `NULL`
#'    in which case they are not used. When `group_colname` and
#'    `id_colname` are defined, then values in `group_colname`
#'    are used to make unique identifiers for each entry in `x`,
#'    and are stored in `id_colname`.
#' @param input_colname `character` string indicating the colname to
#'    use for the input data supplied by `x`. For example when
#'    `input_colname="filename"` then values in `x` are stored in
#'    a column `"filename"`.
#' @param suffix,renameOnes arguments passed to `jamba::makeNames()`,
#'    used when `group_colname` and `id_colname` are defined,
#'    `jamba::makeNames(df[[group_colname]], suffix, renameOnes)`
#'    is used to make unique names for each row.
#' @param colname_hook `function` called on colnames, for example
#'    `jamba::ucfirst()` applies upper-case to the first character
#'    in each colname. When `colname_hook=NULL` then no changes
#'    are made.
#' @param ... additional arguments are passed to `jamba::makeNames()`.
#' 
#' @export
curate_to_df_by_pattern <- function
(x,
 df,
 pattern_colname="pattern",
 group_colname="group",
 id_colname="sample",
 input_colname="filename",
 suffix="_rep",
 renameOnes=TRUE,
 colname_hook=jamba::ucfirst,
 ...)
{
   ## Match pattern with input vector x
   x_match_l <- jamba::provigrep(df[[pattern_colname]], x, returnType="list");
   x_names <- rep(names(x_match_l),
      lengths(x_match_l));
   imatch <- match(x_names,
      df[[pattern_colname]]);
   df_new <- data.frame(check.names=FALSE,
      df[imatch,,drop=FALSE]);
   
   if (length(group_colname) > 0 && 
         length(id_colname) > 0 &&
         group_colname %in% colnames(df)) {
      group_values <- df[imatch, group_colname];
      id_values <- jamba::makeNames(group_values,
         suffix=suffix,
         renameOnes=renameOnes,
         ...);
      df_new[[id_colname]] <- id_values;
      rownames(df_new) <- id_values;
   } else {
      rownames(df_new) <- jamba::makeNames(df_new[[pattern_colname]],
         suffix=suffix,
         renameOnes=renameOnes,
         ...);
   }
   if (length(input_colname) != 1 || any(nchar(input_colname) == 0)) {
      input_colname <- "x";
   }
   df_new[[input_colname]] <- unlist(x_match_l);
   if (length(colname_hook) > 0 && is.function(colname_hook)) {
      colnames(df_new) <- colname_hook(colnames(df_new));
   }
   df_new;
}

