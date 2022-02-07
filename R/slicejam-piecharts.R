
#' Pie chart of feature_type values
#' 
#' Pie chart of feature_type values
#' 
#' This function is called by `slicejam_analysis.Rmd`, and is intended
#' as a simple wrapper to `graphics::pie()` to plot the
#' feature_type values in a `GRanges` object.
#' 
#' @family slicejam plot functions
#' 
#' @param x `character` or `factor` vector that contains feature_type
#'    values to be displayed in the pie chart.
#' @param do_pie `logical` indicating whether to plot the pie chart.
#' @param main `character` string used as the plot title.
#' @param colorSub `character` vector of R colors, whose names are intended
#'    to match values in `x`.
#' @param ... additional arguments are passed to `graphics::pie()`.
#' 
#' @export
feature_to_pie <- function
(x,
 do_pie=TRUE,
 main=NULL,
 colorSub=NULL,
 ...)
{
   j <- table(rmNA(x));
   if (length(j) == 0) {
      return(NULL);
   }
   j_pct <- format(j / sum(j) * 100, digits=1, trim=TRUE);
   j_df <- data.frame(
      j,
      pct=paste0(j_pct, "%")
   );
   main_v <- strsplit(main, "\n")[[1]];
   if (length(main_v) == 4) {
      names(main_v) <- c("Contrast", "Hit", "Signal", "PeakSet");
   } else {
      main_v <- c(Contrast=main, Hit="", Signal="", PeakSet="");
   }
   j_df[,names(main_v)] <- data.frame(as.list(main_v));
   j_names <- names(j);
   colorSub <- rmNA(colorSub[j_names]);
   if (length(colorSub) == length(j_names)) {
      colorSub <- nameVector(rainbowJam(length(j_names)),
         j_names);
   }
   names(j) <- paste0(j_names,
      " (",
      format(trim=TRUE, big.mark=",", j),
      "; ",
      j_pct,
      "%)");
   main <- paste0(main,
      "\n(",
      sum(j),
      " total)");
   if (do_pie) {
      pie(j,
         col=rainbowJam(length(j)));
      title(main=main,
         ...);
   }
   return(j_df);
}

#' Pie chart of feature_type values, wrapper to feature_to_pie()
#' 
#' Pie chart of feature_type values, wrapper to feature_to_pie()
#' 
#' This function is called by `slicejam_analysis.Rmd` and is intended
#' mainly as a wrapper to `feature_to_pie()` in order to prepare
#' multiple pie charts in one multi-panel figure.
#' 
#' Note this function may be replaced with a ggplot2 equivalent, in
#' order to use ggplot2 common color key to help provide more
#' legible labels.
#' 
#' @family slicejam plot functions
#' 
#' @param hit_array `array` output from `se_contrast_stats()` with
#'    three dimensions: hit, contrasts, signal. The arguments
#'    `iHit`, `iCon`, and `iSig` are used to subset the array for
#'    a specific result to use in each pie chart. Each cell of the
#'    array is expected to contain a numeric list with values `c(-1, 1)`
#'    with names that are peak names.
#' @param feature_type_winner `character` or `factor` vector whose names
#'    are expected to be peak names, which are also names of entries
#'    in `hit_array`.
#' @param iHit,iCon,iSig `character` value indicating a specific hit,
#'    contrast, signal, respectively.
#' @param iPeaks `character` vector with specific peaks to use, instead of
#'    using values in `hit_array`.
#' @param main `character` string used as the plot title.
#' @param colorSub `character` vector of R colors, whose names are intended
#'    to match values in `x`.
#' @param ... additional arguments are passed to `feature_to_pie()`.
#' 
#' @export
peak_pie_by_region <- function
(hit_array,
 feature_type_winner=NULL,
 PeakSet="",
 iHit=NULL,
 iSig=NULL,
 iCon=NULL,
 iPeaks=NULL,
 do_pie=TRUE,
 main=NULL,
 colorSub=NULL,
 ...)
{
   ## Wrapper function for peak-to-region
   if (length(iPeaks) > 0) {
      k <- iPeaks;
      if (length(main) == 0) {
         main <- "All Peaks";
      }
      #printDebug("peak_pie_by_region(): ",
      #   "main:", main);
      j_df <- feature_to_pie(feature_type_winner[k],
         main=main,
         colorSub=colorSub,
         do_pie=do_pie);
      return(j_df);
   } else {
      if (length(iHit) == 0) {
         iHit <- dimnames(hit_array)[[1]];
      } else if (is.numeric(iHit)) {
         iHit <- dimnames(hit_array)[[1]][iHit];
      }
      if (length(iSig) == 0) {
         iSig <- dimnames(hit_array)[[3]];
      } else if (is.numeric(iSig)) {
         iSig <- dimnames(hit_array)[[3]][iSig];
      }
      if (length(iCon) == 0) {
         iCon <- dimnames(hit_array)[[2]];
      } else if (is.numeric(iSig)) {
         iCon <- dimnames(hit_array)[[2]][iSig];
      }
      #printDebug("peak_pie_by_region(): ",
      #   "\niCon:", iCon,
      #   "\niHit:", iHit,
      #   "\niSig:", iSig);
      j_df_l <- rbindList(lapply(iHit, function(iHit1) {
         rbindList(lapply(iSig, function(iSig1) {
            rbindList(lapply(iCon, function(iCon1) {
               k <- names(hit_array[iHit1, iCon1, iSig1][[1]]);
               main <- paste(c(iCon1, iHit1, iSig1, PeakSet), collapse="\n");
               j_df <- feature_to_pie(feature_type_winner[k],
                  main=main,
                  colorSub=colorSub,
                  do_pie=do_pie);
               j_df;
            }))
         }))
      }));
      return(j_df_l);
   }
}
