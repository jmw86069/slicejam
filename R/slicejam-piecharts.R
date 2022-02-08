
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
   if (length(main) > 0) {
      main_v <- strsplit(main, "\n")[[1]];
      if (length(main_v) == 4) {
         names(main_v) <- c("Contrast", "Hit", "Signal", "PeakSet");
      } else {
         main_v <- c(Contrast=main, Hit="", Signal="", PeakSet="");
      }
      j_df[,names(main_v)] <- data.frame(as.list(main_v));
      main <- paste0(main,
         "\n(",
         sum(j),
         " total)");
   } else {
      main <- NULL
   }
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


#' ggplot2 Pie chart of peaks by feature_type
#' 
#' ggplot2 Pie chart of peaks by feature_type
#' 
#' @return object of class `c("gg", "ggplot")` that also contains two
#'    attributes:
#'    1. `"panels"`: the unique panel labels used for each facet, delimited
#'    by `"|"`.
#'    2. `"panel_count"`: the number of facet panels, useful when determining
#'    an appropriate figure size. We recommend about 6 inches width and
#'    7 inches height per panel.
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
#' @param iPeaks `character` not used in this implementation.
#' @param do_pie `logical` indicating whether to return the pie chart,
#'    or when `do_pie=FALSE` the `data.frame` is returned.
#' @param colorSub `character` optionally used to define categorical colors.
#' @param label_type `character` indicating the type of pie chart label
#'    to place around the outside: `"text"` uses `ggrepel::geom_text_repel()`,
#'    and `"label"` uses `ggrepel::geom_label_repel()` with colored background.
#' 
#' @export
gg_pie_by_feature_type <- function
(hit_array,
 feature_type_winner,
 PeakSet="",
 iHit=NULL,
 iSig=NULL,
 iCon=NULL,
 iPeaks=NULL,
 add_peak_list=NULL,
 do_pie=TRUE,
 main=NULL,
 colorSub=NULL,
 label_type=c("text", "label"),
 label_size=5,
 base_size=16,
 strip_size=24,
 strip_sub_size=12,
 ...)
{
   #
   label_type <- match.arg(label_type);
   if (length(feature_type_winner) == 0) {
      stop("feature_type_winner must be supplied and non-empty.")
   }
   
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
   
   peak_to_tabledf <- function(k, feature_type_winner, iCon, iHit, iSig)
   {
      j <- table(rmNA(feature_type_winner[k]));
      j_pct <- as.integer(j) / sum(j) * 100;
      j_df <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         feature_type=names(j),
         Count=as.integer(j),
         Scaled=j_pct,
         Percent=paste0(
            format(j_pct,
               digits=1,
               trim=TRUE),
            "%")
      );
      j_df$Offset <- sum(j_df$Count) -(cumsum(j_df$Count) - j_df$Count / 2);
      j_df$Scaled_Offset <- 100 * j_df$Offset / sum(j_df$Count, na.rm=TRUE);
      j_df$Set <- iCon;
      j_df$Hit <- iHit;
      j_df$Signal <- iSig;
      j_df$Label <- paste0(
         j_df$feature_type,
         "\n",
         jamba::formatInt(j_df$Count),
         " peaks\n",
         j_df$Percent);
      j_df;
   }

   jdf <- jamba::rbindList(lapply(iHit, function(iHit1) {
      jamba::rbindList(lapply(iSig, function(iSig1) {
         jamba::rbindList(lapply(iCon, function(iCon1) {
            k <- names(hit_array[iHit1, iCon1, iSig1][[1]]);
            j_df <- peak_to_tabledf(k,
               feature_type_winner,
               iCon1,
               iHit1,
               iSig1);
            j_df;
         }))
      }))
   }));
   #jamba::printDebug("head(jdf, 10):");
   #print(head(jdf, 10));
   
   # optional extra peaks
   if (length(add_peak_list) > 0) {
      jdf_add_list <- lapply(jamba::nameVectorN(add_peak_list), function(iname){
         k <- add_peak_list[[iname]];
         if (is.numeric(k)) {
            k <- names(k)
         }
         j_df <- peak_to_tabledf(k,
            feature_type_winner,
            iCon=iname,
            iHit="",
            iSig="");
      });
      jdf_add <- jamba::rbindList(jdf_add_list);
      #jamba::printDebug("head(jdf_add, 10):");
      #print(head(jdf_add, 10));
      jdf <- jamba::rbindList(list(jdf, jdf_add));
   }
   jdf$Set <- factor(jdf$Set,
      levels=unique(jdf$Set));
   jdf$Hit <- factor(jdf$Hit,
      levels=unique(jdf$Hit));
   jdf$Signal <- factor(jdf$Signal,
      levels=unique(jdf$Signal));
   
   add_span <- function(x, fontsize=strip_sub_size) {
      paste0("<span style='font-size:", fontsize, "pt'>",
         x,
         "</span>")
   }
   jdf$HitLab <- factor(
      add_span(jdf$Hit),
      levels=add_span(levels(jdf$Hit)));
   jdf$SignalLab <- factor(
      add_span(jdf$Signal),
      levels=add_span(levels(jdf$Signal)));
   jdf$Panel <- jamba::pasteByRowOrdered(
      jdf[,c("Set", "Hit", "Signal")],
      sep="|");
   
   if (!do_pie) {
      attr(jdf, "panels") <- levels(jdf$Panel);
      attr(jdf, "panel_count") <- length(levels(jdf$Panel));
      return(jdf);
   }
   
   # Define ggplot2 object
   gg2 <- ggplot2::ggplot(jdf,
      ggplot2::aes(x=1,
         y=Scaled,
         fill=feature_type),
      position="stack") +
      ggplot2::geom_col(color="white",
         position="stack",
         width=1) +
      ggplot2::coord_polar(theta="y",
         direction=-1) +
      ggplot2::scale_x_continuous(
         expand=ggplot2::expansion(mult=c(0, 0.3))) +
      colorjam::theme_jam() +
      ggplot2::theme_void(base_size=base_size) +
      ggplot2::theme(
         #strip.text=ggplot2::element_text(size=strip_size)) +
         strip.text=ggtext::element_markdown(size=strip_size,
            lineheight=1.1)) +
      ggplot2::facet_wrap(~Set + HitLab + SignalLab);
   
   # define categorical colors
   if (length(colorSub) == 0) {
      gg2 <- gg2 +
         colorjam::scale_fill_jam() +
         colorjam::scale_color_jam(invert=TRUE,
            useGrey=7)
   } else {
      gg2 <- gg2 +
         colorjam::scale_fill_manual(values=colorSub) +
         colorjam::scale_color_manual(
            values=jamba::setTextContrastColor(colorSub,
               useGrey=7))
   }
   
   # add text information around the pie chart
   if ("text" %in% label_type) {
      gg2 <- gg2 +
         ggrepel::geom_text_repel(
            segment.color="black",
            color="black",
            ggplot2::aes(x=1.525,
               y=Scaled_Offset,
               label=Label),
            nudge_x=0.9,
            size=label_size,
            show.legend=FALSE)
   } else if ("label" %in% label_type) {
      gg2 +
         ggrepel::geom_label_repel(
            segment.color="black",
            aes(x=1.525,
               color=feature_type,
               y=Scaled_Offset,
               fill=feature_type,
               label=Label),
            nudge_x=0.9,
            size=label_size,
            show.legend=FALSE)
   }
   if (length(main) > 0) {
      gg2 <- gg2 +
         ggplot2::ggtitle(label=main);
   }
   attr(gg2, "panels") <- levels(jdf$Panel);
   attr(gg2, "panel_count") <- length(levels(jdf$Panel));
   gg2;
}
