
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
#' @examples
#' x <- c("NS573_UL3-DexGR_dedup_Single_Fragment_namesort.bam",
#'    "NS573_UL3-VehGR_dedup_Single_Fragment_namesort.bam",
#'    "NS755_UL3-Dex-GR_dedup_SingleFragment_Coordsort.bam",
#'    "NS755_UL3-EtOH-GR_dedup_SingleFragment_Coordsort.bam");
#' 
#' df <- data.frame(
#'    Pattern=c("NS573_UL3-DexGR","NS573_UL3-VehGR","NS755_UL3-Dex-GR","NS755_UL3-EtOH-GR"),
#'    Group=c("UL3_DexGR","UL3_VehGR","UL3_DexGR","UL3_VehGR"),
#'    Batch=c("NS573","NS573","NS755","NS755"),
#'    Label=c("UL3_DexGR_NS573","UL3_VehGR_NS573","UL3_DexGR_NS755","UL3_VehGR_NS755")
#' )
#' df;
#' curate_to_df_by_pattern(x=x, df=df, verbose=TRUE)
#' 
#' curate_to_df_by_pattern(x=x, df=df[,1:3])
#' 
#' @export
curate_to_df_by_pattern <- function
(x,
 df,
 pattern_colname="pattern",
 group_colname="group",
 id_colname=c("label", "sample"),
 input_colname="filename",
 suffix="_rep",
 renameOnes=TRUE,
 colname_hook=jamba::ucfirst,
 sep="_",
 verbose=FALSE,
 ...)
{
   ## Match pattern with input vector x
   pattern_colname <- head(jamba::rmNA(colnames(df)[match(tolower(pattern_colname), 
      tolower(colnames(df)))]), 1);
   if (length(pattern_colname) == 0) {
      pattern_colname <- head(colnames(df), 1);
   }
   group_colname <- jamba::rmNA(colnames(df)[match(tolower(group_colname), 
      tolower(colnames(df)))]);
   id_colname <- head(jamba::rmNA(colnames(df)[match(tolower(id_colname),
      tolower(colnames(df)))]), 1);
   if (length(input_colname) != 1 || any(nchar(input_colname) == 0)) {
      input_colname <- "x";
   }
   if (verbose) {
      jamba::printDebug("curate_to_df_by_pattern(): ",
         "pattern_colname:",
         pattern_colname);
      jamba::printDebug("curate_to_df_by_pattern(): ",
         "group_colname:",
         group_colname);
      jamba::printDebug("curate_to_df_by_pattern(): ",
         "id_colname:",
         id_colname);
      jamba::printDebug("curate_to_df_by_pattern(): ",
         "input_colname:",
         input_colname);
   }
   x_match_l <- jamba::provigrep(df[[pattern_colname]],
      x,
      returnType="list");
   x_names <- rep(names(x_match_l),
      lengths(x_match_l));
   imatch <- match(x_names,
      df[[pattern_colname]]);
   df_new <- data.frame(check.names=FALSE,
      df[imatch,,drop=FALSE]);
   df_new[[input_colname]] <- unlist(x_match_l);
   
   if (length(id_colname) == 0) {
      if (length(group_colname) > 0 &&
            group_colname %in% colnames(df)) {
         label_colnames <- setdiff(colnames(df),
            pattern_colname);
         group_values <- jamba::pasteByRow(df[imatch, label_colnames, drop=FALSE],
            sep=sep,
            ...);
         id_values <- jamba::makeNames(group_values,
            suffix=suffix,
            renameOnes=renameOnes,
            ...);
         id_colname <- "label";
         df_new[[id_colname]] <- id_values;
         rownames(df_new) <- id_values;
      } else {
         id_colname <- input_colname;
         rownames(df_new) <- jamba::makeNames(df_new[[input_colname]],
            suffix=suffix,
            renameOnes=renameOnes,
            ...);
      }
   } else {
      id_values <- jamba::pasteByRow(df[imatch, id_colname, drop=FALSE],
         sep=sep,
         ...);
      rownames(df_new) <- id_values;
   }
   df_colnames <- unique(c(
      setdiff(colnames(df_new), input_colname),
      input_colname));
   df_new <- df_new[, df_colnames, drop=FALSE];
   if (length(colname_hook) > 0 && is.function(colname_hook)) {
      colnames(df_new) <- colname_hook(colnames(df_new));
   }
   df_new;
}

#' Import featureCounts file
#' 
#' @export
import_featurecounts <- function
(file,
 verbose=FALSE,
 rowid_colname="Geneid",
 ...)
{
   if (!file.exists(file)) {
      stop(paste0("File does not exist:",
         file));
   }
   if (length(file) > 1) {
      if (length(names(file)) == 0) {
         names(file) <- jamba::makeNames(file);
      }
      fc_datas <- lapply(file, function(i){
         import_featurecounts(file=i,
            verbose=verbose,
            ...);
      })
   }
   fc_data <- data.table::fread(
      file,
      data.table=FALSE,
      skip=rowid_colname);

   ## Rename first column
   #fc_data[,"Geneid"] <- paste0("peak_",
   #   jamba::padInteger(seq_len((nrow(fc_data)))));
   #rownames(fc_data) <- fc_data[,"Geneid"];
   
   ## Recognize featureCounts columns
   fc_colnames <- intersect(
      c(rowid_colname,
         "Chr",
         "Start",
         "End",
         "Strand",
         "Length"),
      colnames(fc_data));

   ## Recognize counts columns as all other colnames
   filenames <- setdiff(colnames(fc_data),
      fc_colnames);
   attr(fc_data, "fc_colnames") <- fc_colnames;
   attr(fc_data, "filenames") <- filenames;
   return(fc_data);
}
