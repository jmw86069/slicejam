
## volcano plot functions

#' Volcano plot
#' 
#' Draw a volcano plot using a reasonably robust set of default
#' arguments, and with a large number of customization options.
#' The default plot uses smooth scatter plot for much improved
#' display of point density.
#' 
#' This function produces a volcano plot, which consists of
#' change on the x-axis, and significance on the y-axis.
#' 
#' In addition to displaying the volcano plot, this function
#' also displays statistical thresholds, and marks entries
#' as "hits" by up to three conceptual filters:
#' 
#' * "change" - fold change `fold_cutoff`
#' * "significant" - statistical P-value `sig_cutoff`
#' * "detected" - signal `expr_cutoff`
#' 
#' If any cutoff is not defined, that filter is ignored.
#' 
#' Change is usually represented using log2 fold changes,
#' and in this case is labeled using normal scale fold change
#' values. The threshold is defined with `fold_cutoff` using
#' normal space values. The log2 fold change values which
#' have greater magnitude than `fold_cutoff` are marked
#' "changing".
#' 
#' Significance usually represents adjusted P-value, or raw
#' P-value if necessary. The threshold is defined with `sig_cutoff`
#' using a P-value below which entries are marked "significant".
#' 
#' Finally, since some statistical criteria also include a minimum
#' level of signal, a threshold `expr_cutoff` requires an entry to
#' have signal at or above this value to be considered "detected".
#' 
#' The default behavior of `volcano_plot()` is to render a
#' smooth scatter plot. A smooth scatter plot is much more
#' effective at representing the true point density along
#' the figure, which is one of the primary reasons to produce
#' the plot.
#' 
#' ## Highlighting points
#' 
#' The argument `hi_points` can be used to highlight a specific
#' subset of points on the figure, even when `smooth=TRUE`.
#' 
#' Alternatively, `hi_hits=TRUE` will render all statistical
#' hits as points, which will appear on top of the smooth
#' scatter plot when `smooth=TRUE`.
#' 
#' 
#' @param x `data.frame` that contains statistical results with at
#'    least a P-value, and fold change or log2 fold change. It is
#'    useful to contain a column with mean expression, and a column
#'    with a relevant label.
#' @param n `integer` indicating the number of subset points to plot
#'    for testing purposes.
#' @param lfc_colname `character` string or vector used to match
#'    `colnames(x)` whose values should be log2 fold changes.
#'    A direct match to `colnames(x)` is performed
#'    first, then if no column is found, the values are used as
#'    regular expression patterns in order until the first
#'    matching colname is found. Note that `lfc_colname` is used
#'    in preference to `fold_colname`.
#'    The colname used will appear as the x-axis label.
#' @param fold_colname `character` string or vector used to match
#'    `colnames(x)` whose values should be fold changes. Note that
#'    if `lfc_colname` successfully finds a value, the `fold_colname`
#'    is not used.
#'    The colname if used will appear as the x-axis label.
#' @param sig_colname `character` string or vector used to match
#'    `colnames(x)` whose values should contain P-values of significance.
#'    The P-values can be unadjusted (raw) P-values, or adjusted
#'    P-values. The P-values are expected not to be `-log10()`
#'    transformed.
#'    The colname used will appear as the y-axis label.
#' @param expr_colname `character` string or vector used to match
#'    `colnames(x)` whose values should contain expression mean values.
#'    This column is only used when `expr_cutoff` is defined and
#'    is applied to the filter criteria for statistical hits.
#' @param label_colname `character` string or vector used to match
#'    `colnames(x)` whose values should contain a useful label,
#'    for example gene symbol or assay identifier.
#' @param sig_cutoff `numeric` threshold for values in `sig_colname`,
#'    where values at or below `sig_cutoff` can be considered
#'    statistically significant.
#' @param fold_cutoff `numeric` threshold for values in `lfc_colname`
#'    or `fold_cutoff`, where normal fold change values at or above
#'    `fold_cutoff` can be considered statistically significant.
#'    Note that when `lfc_colname` is being used, its values are
#'    converted to normal fold change before applying this filter.
#' @param expr_cutoff `numeric` threshold for values in `expr_colname`
#'    when `expr_colname` is defined, where values in `expr_colname`
#'    at or above `expr_cutoff` can be considered statistically
#'    significant. This threshold is useful to filter out potential
#'    statistical hits whose signal is below a noise signal threshold.
#' @param main `character` string used as the main title of the figure.
#' @param submain `character` string used as a sub-title of the figure.
#' @param symmetric_axes `logical` indicating whether the x-axis
#'    log fold change range should be symmetric above and below zero.
#' @param do_cutoff_caption `logical` indicating whether to include
#'    a text caption with the statistical cutoff values used.
#' @param caption_cex `numeric` caption font size adjustment.
#' @param sig_floor `numeric` value indicating the lowest P-value to
#'    display on the y-axis, where values below `sig_floor` are
#'    fixed to the top of the y-axis range. This value is useful
#'    when extremely low P-values (1e-100) compress the useful
#'    range of P-values on the y-axis.
#' @param min_sig_range `numeric` value indicating a y-axis P-value
#'    the must be included in the visual range. This argument is useful
#'    when P-values are not very significant, and you want to make
#'    sure the y-axis range shows a minimum amount of the significant
#'    region to be visually interpretable in that context.
#' @param fold_ceiling `numeric` value indicating the maximum fold
#'    change range displayed on the x-axis, useful to prevent extremely
#'    large fold changes from compressing the useful visible range
#'    of the figure.
#' @param min_fold_range `numeric` indicating the lowest fold change
#'    to include in the x-axis visual range. This argument is useful
#'    when fold changes are low and the x-axis range would otherwise
#'    be too small to be very useful.
#' @param n_x_labels,n_y_labels `integer` used by `pretty()` to determine
#'    the approximate number of x-axis and y-axis labels to display,
#'    respectively.
#' @param xlim,ylim `numeric` used to define specific `xlim` and `ylim`
#'    axis ranges. When `NULL` the ranges are defined automatically,
#'    and when not `NULL` the `xlim` or `ylim` values are used directly.
#' @param pt_cex,pt_pch `numeric` used to define point size and shape,
#'    used only when individual points are displayed.
#' @param hit_type `character` string used to label points that meet
#'    the statistical cutoffs as `"hits"`, but where it may be useful
#'    to indicate the type of entry being tested. For example,
#'    `hit_type="genes"` indicates that each row represents a gene;
#'    `hit_type="probes"` indicates each row represents a probe;
#'    `hit_type="transcripts"` indicates each row represents a transcript.
#' @param color_set `character` vector of R colors, used only when individual
#'    points are display. The names override default values, and may include:
#'    * `"base"` - the base color of all points on the plot
#'    * `"up"` - the color for up-regulated points that meet all
#'    statistical cutoffs to be a "hit".
#'    * `"down"` - the color for down-regulated points that meet all cutoffs
#'    * `"hi"` - base color for highlighted points, used when `hi_points`
#'    is defined.
#'    * `"hi_up"` - color for highlighted up-regulated points.
#'    * `"hi_down"` - color for highlighted down-regulated points.
#' @param border_set `NULL` or `character` vector of R colors, used to
#'    define point border colors such as `pch=21` which is a filled circle
#'    with border. When `border_set=NULL` then it is defined by
#'    `jamba::makeColorDarker(color_set)`.
#' @param point_colors,border_colors optional `character` vector of R colors
#'    recycled to length `nrow(x)`, used to specify the exact color of each
#'    point in `x`. This argument is useful to colorize certain specific
#'    points that may otherwise not meet statistical criteria.
#' @param abline_color `character` string with R color used to color
#'    the abline that indicates the x-axis `fold_cutoff` value, and
#'    y-axis `sig_cutoff` value.
#' @param smooth `logical` indicating whether points should be drawn
#'    as a smooth scatter plot, using `jamba::plotSmoothScatter()`.
#'    When `smooth=FALSE` individual points are drawn, using
#'    `point_colors`, or when `point_colors` is not defined the
#'    default is to use `color_set` to colorize points based upon
#'    statistical cutoffs.
#' @param smooth_func `function` used to plot points when `smooth=TRUE`,
#'    by default `jamba::plotSmoothScatter()` which has some benefits
#'    over default `graphics::smoothScatter()`.
#' @param smooth_ramp `character` vector of R colors which defines
#'    the color gradient to use when `smooth=TRUE`.
#' @param blockarrow `logical` indicating whether block arrows should
#'    be displayed and used to indicate the number of statistical hits.
#' @param blockarrow_colors,blockarrow_font,blockarrow_label,blockarrow_shadowtext
#'    arguments used when `blockarrow=TRUE`.
#' @param tophist `logical` indicating whether to display a histogram
#'    at the top of the volcano plot figure.
#' @param tophist_cutoffs,tophist_breaks,tophist_fraction,tophist_by
#'    arguments used when `tophist=TRUE`.
#' @param hi_points `character` vector indicating points to highlight
#'    in the volcano plot, where values should match `rownames(x)`.
#'    This argument is useful to highlight a specific subset of points of
#'    interest on the figure. Note that `hi_points` are always
#'    rendered as individual points even when `smooth=TRUE`.
#' @param hi_hits `logical` indicating whether rows that meet all
#'    statistical cutoffs and are considered "hits" should also be
#'    treated as `hi_points` for the purpose of rendering individual
#'    points.
#' @param hi_cex `numeric` size adjustment for highlight points,
#'    relative to the size of other points in the figure.
#' @param do_both `logical` indicating whether to draw both a smooth
#'    scatter and individual points on the same figure.
#' @param label_hits `logical` indicating whether to add a text label
#'    for points that are statistical hits.
#' @param add_plot `logical` indicating whether the plot should be
#'    added to an existing plot, or when `add_plot=FALSE` a new
#'    plot is created. This argument is useful to re-run the
#'    same volcano plot with alternate parameters, for example
#'    to display different subsets of highlighted points.
#' @param xlab,ylab `character` strings used to specify the exact
#'    x-axis label and y-axis label. When either value is `NULL`
#'    the default is to use the relevant colname: x-axis uses
#'    either `lfc_colname` or `fold_colname`; y-axis uses `sig_colname`.
#' @param cex.axis `numeric` adjustment for axis label font sizes.
#' @param include_axis_prefix `logical` indicating whether to include
#'    a prefix for the x-axis and y-axis labels: x-axis `"Change"`;
#'    y-axis `"Significance"`.
#' @param `numeric` vector used to ensure that each margin size is
#'    at least a minimum value, applied to `par("mar")` via
#'    the function `pmax()`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param transformation `function` passed to `smooth_func` used to
#'    adjust the visual contrast of the resulting density plot.
#' @param nbin `numeric` value passed to `smooth_func` and used
#'    by `jamba::plotSmoothScatter()` to adjust the number of
#'    bins used to display the density of points, where a higher
#'    value shows more detail, and a lower value shows less detail.
#'    
#' @examples
#' n <- 150;
#' set.seed(12);
#' x_lfc <- (rnorm(n) * 1);
#' x_lfc <- x_lfc^2 * sign(x_lfc);
#' x_lfc <- x_lfc[order(-abs(x_lfc) + rnorm(n) / 2)];
#' x_pv <- sort(10^-(rnorm(n)*1.5)^2);
#' x <- data.frame(
#'    Gene=paste("gene", seq_len(n)),
#'    `log2fold Group-Control`=x_lfc,
#'    `P.Value Group-Control`=x_pv[order(-abs(x_lfc))],
#'    check.names=FALSE);
#' volcano_plot(x);
#' 
#' par("mfrow"=c(2, 1));
#' volcano_plot(x);
#' volcano_plot(x);
#' par("mfrow"=c(1, 1));
#' 
#' x[["fold Group-Control"]] <- log2fold_to_fold(x[["log2fold Group-Control"]]);
#' x[["adj.P.Val Group-Control"]] <- x[["P.Value Group-Control"]];
#' 
#' volcano_plot(x, hi_hits=TRUE);
#' 
#' @export
volcano_plot <- function
(x,
 n=NULL,
 lfc_colname=c("logfc", "log2fold", "log2fc", "lfc", "logratio", "log2ratio"),
 fold_colname=c("fold", "fc", "ratio"),
 sig_colname=c("adj.P.Val", "P.Value"),
 expr_colname=c("mgm", "groupmean", "mean", "AveExpr", "fkpm", "tpm"),
 label_colname=c("gene", "probe"),
 sig_cutoff=0.05,
 fold_cutoff=1.5,
 expr_cutoff=NULL,
 main="Volcano Plot",
 submain=NULL,
 symmetric_axes=TRUE,
 do_cutoff_caption=TRUE,
 caption_cex=0.8,
 sig_floor=1e-20,
 min_sig_range=1e-4,
 fold_ceiling=64,
 min_fold_range=8,
 n_x_labels=12,
 n_y_labels=7,
 xlim=NULL,
 ylim=NULL,
 pt_cex=0.9,
 pt_pch=21,
 hit_type="hits",
 color_set=c(base="#77777777",
    up="#99000088",
    down="#00009988",
    hi="#FFDD55FF",
    hi_up="#FFDD55FF",
    hi_down="#FFDD55FF"),
 border_set=NULL,
 point_colors=NULL,
 border_colors=NULL,
 abline_color="#000000AA",
 smooth=TRUE,
 smooth_func=jamba::plotSmoothScatter,
 smooth_ramp=colorRampPalette(c("white", "lightblue", "lightskyblue3", "royalblue", "darkblue", "orange", "darkorange1", "orangered2")),
 blockarrow=TRUE,
 blockarrow_colors=c(hit="#E67739FF",
    up="#990000FF",
    down="#000099FF"),
 blockarrow_font=1,
 blockarrow_cex=c(1.2, 1.2),
 blockarrow_label_cex=1,
 #blockarrow_label_color="#FFFFAA", blockArrowLabelUpColor="#FFFFAA", blockArrowLabelDownColor="#FFFFAA",
 blockarrow_shadowtext=TRUE,
 tophist=FALSE,
 tophist_cutoffs=c("pvalue", "foldchange"),
 tophist_breaks=100,
 tophist_color="#000099FF",
 tophist_fraction=1/3,
 tophist_by=0.20,
 hi_points=NULL,
 hi_hits=FALSE,
 hi_cex=1,
 do_both=FALSE,
 label_hits=FALSE,
 add_plot=FALSE,
 xlab=NULL,
 ylab=NULL,
 cex.axis=1.2,
 include_axis_prefix=TRUE,
 mar_min=c(5, 4, 6, 5),
 verbose=FALSE,
 transformation=function(x){x^0.24;}, # use 0.14 for very large datasets
 nbin=256,
 ...)
{

   ## Purpose is to wrapper a simple volcano plot with log-friendly axis labels
   ##
   ## orientation can be "up" for upright, "right" to tilt 90 degree on its side
   ##
   ## minPvalue is designed to set NA or '0' values to a fixed value
   ## pvalueFloor is designed to set extremely low P-values (1e-150) to a minimum
   ## sectionLabelSpacing is the fraction of the y-axis range away from each border
   ##    used to position the highlight hit counts labels (if highlightPoints is not NULL)
   ##
   ## A change was made to use pointColorSet to set the colors, instead of using
   ## hitColorSet, so the colors could independently be manipulated.
   ## To use this scheme, set usePointColorSet=TRUE
   ##
   ## The graph plots fold change on the x-axis, but with log2 scale, the numbers
   ## are technically being reported in normal space. However, to give it a log2 axis label:
   ## xlab=expression(log[2]*"(Fold change)"
   blockarrow_cex <- rep(blockarrow_cex, length.out=2);
   
   ## Only alter parameters which were provided in the function call,
   ## but keep default values which were not defined
   color_set <- update_function_params(function_name="volcano_plot",
      param_name="color_set",
      new_values=color_set);
   border_set_default <- jamba::makeColorDarker(color_set,
      darkFactor=2,
      sFactor=1);
   border_set_default <- jamba::alpha2col(border_set_default,
      alpha=pmin(jamba::col2alpha(border_set_default) + 0.5, 1))
   border_set <- update_list_elements(source_list=border_set_default,
      update_list=border_set,
      verbose=FALSE);
   if (verbose) {
      jamba::printDebug("volcano_plot(): ",
         "color_set:");
      jamba::printDebugI(color_set);
      jamba::printDebug("volcano_plot(): ",
         "border_set:");
      jamba::printDebugI(border_set);
   }

   ## Allow processing only a subset 'n' rows of data
   if (length(n) == 0) {
      n <- nrow(x);
   } else if (n < nrow(x)) {
      x <- head(x, n);
   }
   
   if (length(rownames(x)) == 0 || length(jamba::tcount(rownames(x), minCount=2)) > 0) {
      rownames(x) <- jamba::makeNames(rep("row", n));
   }
   
   ## label
   label_colname <- find_colname(label_colname,
      x);
   
   ## significance
   sig_colname <- find_colname(sig_colname,
      x,
      col_types=c("numeric", "integer"));
   if (length(sig_colname) == 0) {
      stop("sig_colname was not defined.");
   }
   sig_values <- x[[sig_colname]];
   ## consider replacing NA with 1
   #sig_values[is.na(sig_values)] <- 1;
   if (verbose) {
      jamba::printDebug("volcano_plot(): ",
         "sig_colname:",
         sig_colname);
   }
   met_sig <- rep(TRUE, nrow(x));
   if (length(sig_cutoff) > 0) {
      met_sig <- (!is.na(sig_values) &
            sig_values <= head(sig_cutoff, 1));
   }
   
   ## fold change
   lfc_colname <- find_colname(lfc_colname,
      x,
      col_types=c("numeric", "integer"));
   fold_colname <- find_colname(fold_colname,
      x,
      col_types=c("numeric", "integer"),
      exclude_pattern=c("log", "lfc"));
   met_fold <- rep(TRUE, nrow(x));
   if (length(lfc_colname) > 0) {
      lfc_values <- x[[lfc_colname]];
      fold_colname_label <- lfc_colname;
   } else if (length(fold_colname) > 0) {
      lfc_values <- fold_to_log2fold(x[[fold_colname]]);
      fold_colname_label <- fold_colname;
   } else {
      stop("Must provide either fold_colname or lfc_colname.");
   }
   ## consider replacing NA lfc_values with 0?
   #lfc_values[is.na(lfc_values)] <- 0;
   if (length(fold_cutoff) > 0 && head(fold_cutoff, 1) > 1) {
      met_fold <- (!is.na(x[[lfc_colname]]) &
            abs(x[[lfc_colname]]) >= log2(head(fold_cutoff, 1)));
   } else {
      fold_cutoff <- NULL;
   }

   ## expression
   met_expr <- rep(TRUE, nrow(x));
   if (length(expr_cutoff) > 0 && length(expr_colname) > 0) {
      expr_cutoff <- head(expr_cutoff, 1);
      expr_colname <- find_colname(expr_colname,
         x,
         max=Inf,
         col_types=c("integer", "numeric"));
      if (length(expr_colname) > 0) {
         if (verbose) {
            jamba::printDebug("volcano_plot(): ",
               "applying expr_cutoff:",
               format(expr_cutoff, digits=3),
               " to expr_colname:",
               expr_colname);
         }
         expr_max <- matrixStats::rowMaxs(
            as.matrix(x[, expr_colname, drop=FALSE]),
            na.rm=TRUE);
         met_expr <- (!is.na(expr_max) &
               expr_max >= expr_cutoff);
         if (verbose) {
            jamba::printDebug("volcano_plot(): ",
               sum(met_expr),
               " values of ",
               length(met_expr),
               " met the threshold.");
         }
      }
   } else {
      expr_cutoff <- NULL;
   }
   
   hits_both <- (met_sig & met_fold & met_expr);
   hits_up <- (lfc_values > 0 & hits_both);
   hits_dn <- (lfc_values < 0 & hits_both);
   pv_hits_up <- sum(hits_up);
   pv_hits_dn <- sum(hits_dn);
   pv_hits <- sum(hits_both);
   
   #pvHitsUpWhich <- which(x[,pvCol] <= pvalueCutoff & (x[,fcCol]) >= fcCutoff & metIntensity);
   #pvHitsDownWhich <- which(x[,pvCol] <= pvalueCutoff & -(x[,fcCol]) >= fcCutoff & metIntensity);
   #pvHitsBothWhich <- which(x[,pvCol] <= pvalueCutoff & abs(x[,fcCol]) >= fcCutoff & metIntensity);
   #pvHitsUp <- length(pvHitsUpWhich);
   #pvHitsDown <- length(pvHitsDownWhich);
   #pvHits <- length(pvHitsBothWhich);
   #upHits <- rownames(x)[pvHitsUpWhich];
   #downHits <- rownames(x)[pvHitsDownWhich];
   
   ## Allow highlighting a subset of points
   if (length(hi_points) > 0) {
      if (is.numeric(hi_points)) {
         hi_points <- seq_len(nrow(x)) %in% hi_points;
      } else {
         hi_points <- rownames(x) %in% hi_points;
      }
   } else {
      hi_points <- rep(FALSE, nrow(x));
   }
   if (hi_hits) {
      hi_points <- hits_both | hi_points;
   }
   if (any(!is.na(hi_points) & hi_points) && length(label_colname) > 0) {
      hi_labels <- x[[label_colname]][hi_points];
   }
   
   ## Define point colors by point_type
   point_type <- rep("base", nrow(x));
   point_type[!hits_both & hi_points] <- "hi";
   point_type[hits_up & !hi_points] <- "up";
   point_type[hits_up & hi_points] <- "hi_up";
   point_type[hits_dn & !hi_points] <- "down";
   point_type[hits_dn & hi_points] <- "hi_down";
   if (verbose) {
      print(table(point_type));
   }
   if (length(point_colors) == 0) {
      point_colors <- color_set[point_type];
   } else {
      point_colors <- rep(point_colors, length.out=nrow(x));
   }
   if (length(border_colors) == 0) {
      border_colors <- border_set[point_type];
   } else {
      border_colors <- rep(border_colors, length.out=nrow(x));
   }
   if (verbose) {
      jamba::printDebug("volcano_plot(): ",
         "head(unique(point_type)):");
      jamba::printDebugI(head(unique(point_type), 20));
      jamba::printDebug("volcano_plot(): ",
         "head(unique(point_colors)):");
      jamba::printDebugI(head(unique(point_colors), 20));
      jamba::printDebug("volcano_plot(): ",
         "head(unique(border_colors)):");
      jamba::printDebugI(head(unique(border_colors), 20));
   }
   
   ## y-axis values
   if (length(sig_floor) == 1 &&
         sig_floor > 0 &&
         any(!is.na(sig_values) & sig_values < sig_floor)) {
      sig_values[sig_values < sig_floor] <- sig_floor;
   }
   y_values <- -log10(sig_values);
   
   ## x-axis values
   if (length(fold_ceiling) == 1 && !is.na(fold_ceiling)) {
      lfc_cap <- (!is.na(lfc_values) & abs(lfc_values) > log2(fold_ceiling));
      if (any(lfc_cap)) {
         lfc_values[lfc_cap] <- log2(fold_ceiling) * sign(lfc_values[lfc_cap]);
      }
   }
   x_values <- lfc_values;

   #if (!is.null(pointsToLabel)) {
   #   row2name <- c(jamba::nameVector(rownames(x), x[,geneColumn]),
   #      jamba::nameVector(rownames(x)));
   #   pointsToLabel <- jamba::flipVector(row2name[pointsToLabel],
   #      makeNamesFunc=c);
   #}
   if (any(!is.na(hi_points) & hi_points)) {
      if (verbose) {
         jamba::printDebug("Creating subset of points for highlighting.");
      }

      ## Generate hit counts among the highlighted points
      hi_up_count <- sum(point_type[hi_points] %in% "hi_up");
      hi_dn_count <- sum(point_type[hi_points] %in% "hi_dn");
      hi_other_count <- sum(hi_points) - hi_up_count - hi_dn_count;
   }
   ## optional plot_only_subset=TRUE here

   ## optional do_jitter here

   if (length(xlim) == 0) {
      fold_min <- log2(max(c(fold_cutoff, 2)) * 1.3) * c(-1, 1);
      min_fold_range <- log2(head(abs(min_fold_range), 1)) * c(-1, 1);
      xlim <- range(c(x_values,
         min_fold_range,
         fold_min),
         na.rm=TRUE);
   }
   if (symmetric_axes) {
      xlim <- max(abs(xlim)) * c(-1, 1);
   }
   if (length(sig_floor) == 0 || any(!is.na(sig_floor) & sig_floor < 1e-300)) {
      sig_floor <- 1e-300;
   }
   if (length(ylim) == 0) {
      sig_min <- max(-log10(c(sig_cutoff, 0.05))) * 1.3;
      # min_sig_range
      ylim <- range(c(0,
         min(c(max(y_values),
            -log10(sig_floor)), na.rm=TRUE),
         -log10(min_sig_range)),
         na.rm=TRUE);
   }

   if (any(!is.na(hi_points) & hi_points)) {
      ###############################################
      ## Position labels in quadrants of the plot
      namesX <- c(xlim[1] - xlim[1]/12,
         xlim[2] - xlim[2]/12);
      namesY <- rep(ylim[2] - ylim[2]/12,
         2);
      namesLabel <- paste(
         c(format(big.mark=",", hi_dn_count),
            format(big.mark=",", hi_up_count)),
         hit_type);
      namesLabel <- gsub("^(1 .+)s$", "\\1", namesLabel);
      if (verbose) {
         jamba::printDebug("volcano_plot(): ",
            "namesLabel: ",
            namesLabel);
      }
   }
   
   ## Overall Title
   do_overall_title <- function() {
      if ((length(main) > 0 && nchar(main) > 0) ||
            (length(submain) > 0 && nchar(submain) > 0)) {
         origPar1 <- par("xpd"=TRUE);
         font.main <- 1;
         if (length(main) > 0 && nchar(main) > 0) {
            nlines_main <- lengths(strsplit(main, "\n"));
            line_main <- parMar[3] - 0.5 - 1.2 * nlines_main;
            title(main=main,
               line=line_main,
               font.main=font.main,
               cex.main=1.5);
         }
         if (length(submain) > 0 && nchar(submain) > 0) {
            title(main=submain,
               line=line_main - nlines_main / 2 - 0.4,
               cex.main=1,
               font.main=font.main);
         }
         par(origPar1);
      }
   }
   
   parList <- list();
   #parList[["prePlots"]] <- origPar;

   
   ## labelCoords will have the return data from addNonOverlappingLabels() but only
   ## if we end up calling that method
   labelCoords <- NULL;
   parMarXpd <- par("mar", "xpd");
   parMar <- parMarXpd$mar;
   jamba::printDebug("parMar: ", parMar);
   on.exit(par(parMarXpd));
   if (length(mar_min) > 0) {
      parMar <- pmax(parMar, mar_min);
   }

   ##########################################################
   ## tophist - histogram of hits along top border
   if (!add_plot && length(tophist) > 0 && is.logical(tophist) && tophist) {
      if (length(grep("pval", tophist_cutoffs)) > 0) {
         ## Apply P-value filtering
         pvHitsWhich <- met_sig;
      } else {
         pvHitsWhich <- rep(TRUE, nrow(x));
      }
      if (length(grep("fc|fold", tophist_cutoffs)) > 0) {
         ## Apply fold change filtering
         fcHitsWhich <- met_fold;
      } else {
         fcHitsWhich <- 1:n;
      }
      hist_which <- (pvHitsWhich & fcHitsWhich);

      ## Color the bars consistent with the block arrows
      tophist_colors <- blockarrow_colors[c("down", "up")];
      tophist_breaks1 <- rev(seq(from=-log2(fold_cutoff), to=xlim[1], by=-tophist_by));
      tophist_breaks2 <- seq(from=log2(fold_cutoff), to=xlim[2], by=tophist_by);
      if (!xlim[1] %in% tophist_breaks1) {
         tophist_breaks1 <- c(tophist_breaks1[1] - tophist_by,
            tophist_breaks1);
      }
      if (!xlim[2] %in% tophist_breaks2) {
         tophist_breaks2 <- c(tophist_breaks2,
            tail(tophist_breaks1, 1) + tophist_by);
      }
      tophist_breaks <- unique(c(tophist_breaks1, tophist_breaks2));
      tophist_data <- hist(x_values[hist_which],
         breaks=tophist_breaks,
         plot=FALSE);
      plot_zones <- matrix(c(1,2),
         nrow=2);
      layout(plot_zones,
         widths=c(1),
         heights=c(tophist_fraction,
            1-tophist_fraction));
      ## Change margins, then plot the top histogram
      ## default is par(mar=c(5.1, 4.1, 4.1, 2.1));
      par("mar"=c(1, parMar[2], parMar[3], parMar[4]));

      #parList[["preTopHist"]] <- origPar;
      #parList[["preTopHist"]]$mar <- c(1, parMar[2], parMar[3], parMar[4]);
      
      r1 <- as.integer(length(tophist_breaks)/2);
      tophist_col <- rep(tophist_colors, c(r1, r1+1));
      parAxs <- par("xaxs"="i", "yaxs"="i");
      expand <- c(0.04, 0.04);
      xlim4 <- sort((c(-1,1) * diff(xlim) * expand[1]/2) + xlim);
      graphics:::plot.histogram(tophist_data,
         freq=TRUE,
         xlim=xlim4,
         ylim=range(tophist_data$counts) + c(0, 1),
         main="", xlab="", ylab="",
         axes=FALSE,
         col=tophist_col);
      jamba::minorLogTicksAxis(1,
         logBase=2,
         displayBase=2,
         majorCex=cex.axis,
         minorCex=cex.axis*0.7,
         symmetricZero=TRUE,
         doLabels=FALSE,
         doMinorLabels=FALSE,
         offset=0,
         ...);
      box();
      par(parAxs);
      prettyAt1 <- pretty(c(0, max(tophist_data$counts) + 1), n=5);
      prettyAt1 <- prettyAt1[prettyAt1 <= max(tophist_data$counts)];
      axis(2,
         las=2,
         cex.axis=1.3,
         at=prettyAt1,
         ...);
      title(ylab="Hits");
      #parList[["postTopHist"]] <- origPar;
      #parList[["postTopHist"]]$mar <- c(parMar[1], parMar[2], 2, parMar[4]);
      do_overall_title();
      par("mar"=c(parMar[1], parMar[2], 3, parMar[4]));
   } else {
      par("mar"=parMar);
   }

   if (length(xlab) == 0) {
      if (include_axis_prefix) {
         xlab <- paste0("Change (", fold_colname_label, ")");
      } else {
         xlab <- fold_colname_label;
      }
   }
   if (length(ylab) == 0) {
      if (include_axis_prefix) {
         ylab <- paste0("Significance (", sig_colname, ")");
      } else {
         ylab <- sig_colname;
      }
   }
   
   ######################
   ## Smooth scatter plot
   if (do_both) {
      smooth <- TRUE;
   }
   if (smooth) {
      ## Upright volcano plot (standard orientation)
      smooth_func(x=x_values,
         y=y_values,
         xaxt="n",
         yaxt="n",
         xlab="",
         ylab="",
         transformation=transformation,
         use_raster=TRUE,
         xlim=xlim,
         ylim=ylim,
         nbin=nbin,
         colramp=smooth_ramp,
         add=FALSE,
         nrpoints=0,
         ...);
      title(xlab=xlab,
         line=mean(c(2.5, parMar[1] - 2.5)),
         cex.lab=cex.axis,
         ...);
      title(ylab=ylab,
         line=mean(c(3, parMar[2] - 1.2)),
         cex.lab=cex.axis,
         ...);
      parUsr <- par("usr");
      #parList[["postSmoothScatter"]] <- par(no.readonly=TRUE);
   }
   
   
   ## Non-smooth scatter points
   if (!smooth) {
      ## quick blank plot to set axis ranges
      plot(NULL,
         xlim=xlim,
         ylim=ylim,
         xaxt="n",
         yaxt="n",
         xlab="",
         ylab="",
         ...);
      title(xlab=xlab,
         line=mean(c(2.5, parMar[1] - 2.5)),
         #line=2.5,
         cex.axis=cex.axis,
         ...);
      title(ylab=ylab,
         line=mean(c(3, parMar[2] - 1.2)),
         #line=4,
         cex.axis=cex.axis,
         ...);
   }
   if (any(!is.na(hi_points) & hi_points)) {
      if (do_both) {
         points(x=x_values[!hi_points],
            y=y_values[!hi_points],
            pch=pt_pch,
            cex=pt_cex,
            bg=point_colors[!hi_points],
            col=border_colors[!hi_points],
            xaxt="n",
            yaxt="n",
            ...);
      }
      points(x=x_values[hi_points],
         y=y_values[hi_points],
         pch=pt_pch,
         cex=hi_cex,
         bg=point_colors[hi_points],
         col=border_colors[hi_points],
         xaxt="n",
         yaxt="n",
         ...);
      ## Optionally labels the highlighted points, using the addNonOverlappingLabels() function
      if (label_hits) {
         labelCoords <- addNonOverlappingLabels(
            x=x_values[hi_points],
            y=y_values[hi_points],
            initialAngle=initialAngle,
            txt=hi_labels,
            n=labelN,
            labelCex=labelCex,
            labelMar=labelMar,
            boxColor=boxColor,
            boxBorderColor=boxBorderColor,
            initialRadius=initialRadius,
            fixedCoords=labelFixedCoords,
            ...);
      }
      ## Add text labels indicating the number of highlighted points
      textAdjX <- c(0, 1);
      if (label_hits) {
         #jamba::drawLabels(preset="topleft",
         for(i in seq_along(namesLabel)) {
            text(x=namesX[i],
               y=namesY[i],
               label=namesLabel[i],
               adj=c(textAdjX[i], 0.5));
         }
      }
   } else if (!smooth || do_both) {
      if (verbose) {
         jamba::printDebug("volcano_plot(): ",
            "Following do_both=TRUE without hi_points.");
      }
      points(x=x_values,
         y=y_values,
         pch=pt_pch,
         cex=pt_cex,
         bg=point_colors,
         col=border_colors,
         xaxt="n",
         yaxt="n",
         ...);
   }
   #parList[["postScatter"]] <- par(no.readonly=TRUE);

   multiGenesUp <- character(0);
   multiGenesDown <- character(0);
   if (!add_plot) {
      ## Define labels
      label_up <- paste(format(big.mark=",", sum(hits_up)),
         hit_type);
      if (sum(hits_up) == 1) {
         label_up <- gsub("s$", "", label_up);
      }
      label_up <- paste(label_up, "up");
      
      label_dn <- paste(format(big.mark=",", sum(hits_dn)),
         hit_type);
      if (sum(hits_dn) == 1) {
         label_dn <- gsub("s$", "", label_dn);
      }
      label_dn <- paste(label_dn, "down");
      
      label_both <- paste0(format(big.mark=",", sum(hits_both)),
         " significant ",
         hit_type);
      if (sum(hits_both) == 1) {
         label_both <- gsub("s$", "", label_both);
      }
      
      ## optional hits_by_gene

      ## draw hit labels   
      if (!add_plot && label_hits && !blockarrow) {
         jamba::drawLabels(preset=c("topleft", "topright"),
            txt=c(label_dn, label_up),
            labelCex=1,
            drawBox=FALSE);
      }
      
      ## Block Arrows
      #parList[["preBlockArrows"]] <- par(no.readonly=TRUE);
      #origPar1 <- par("xpd"=FALSE);
      #on.exit(par(origPar1), add=TRUE);
      
      y_at <- unique(as.integer(pretty(ylim, n=n_y_labels)));
      logAxis(2, 
         at=unique(as.integer(pretty(ylim, n=n_y_labels))),
         value=FALSE, 
         base=10, 
         makeNegative=TRUE, 
         cex.axis=cex.axis*1,
         ...);
      ## Add small label indicating the threshold
      if (!-log10(sig_cutoff) %in% y_at) {
         axis(2,
            at=-log10(sig_cutoff),
            labels=format(sig_cutoff),
            las=2,
            cex=cex.axis*0.8)
      }
      #logAxis(1, at=unique(as.integer(pretty(xRange, n=nXlabels))),
      #   value=TRUE, base=2, cex.axis=cex.axis*0.8, ...);
      jamba::minorLogTicksAxis(1,
         logBase=2,
         displayBase=2,
         majorCex=cex.axis,
         minorCex=cex.axis*0.7,
         symmetricZero=TRUE,
         offset=0,
         ...);
      par("xpd"=FALSE);
      if (length(sig_cutoff) > 0) {
         abline(h=-log10(sig_cutoff),
            lty="dashed", 
            col=abline_color);
      }
      if (length(fold_cutoff) > 0) {
         abline(v=unique(log2(fold_cutoff) * c(-1, 1)),
         lty="dashed", 
         col=abline_color);
      }
      if (blockarrow) {
         par("xpd"=TRUE);
         hitCol <- hsv(h=0.06, 
            s=0.75, 
            v=0.9, 
            alpha=1);
         if (tophist) {
            right_adj <- 1.2;
         } else {
            right_adj <- 1;
         }
         blockarrow_label_colors <- jamba::setTextContrastColor(
            jamba::makeColorDarker(blockarrow_colors,
               darkFactor=1.3,
               sFactor=1.3))

         blockarrow_cex <- rep(blockarrow_cex, length.out=2);
         blockArrowMargin(axisPosition="rightAxis",
            ybottom=-log10(sig_cutoff),
            arrowPosition="top",
            doShadowText=blockarrow_shadowtext,
            blockWidthPercent=5*blockarrow_cex[2]*right_adj,
            arrowLabel=label_both,
            col=blockarrow_colors[["hit"]],
            labelCex=blockarrow_label_cex * blockarrow_cex[2],
            labelFont=blockarrow_font, 
            arrowLabelColor=blockarrow_label_colors[["hit"]],
            arrowLabelBorder=jamba::alpha2col(blockarrow_label_colors[["hit"]], 0.3));
         if (length(fold_cutoff) == 0) {
            fold_cutoff <- 1;
         }
         blockArrowMargin(axisPosition="topAxis",
            xleft=log2(fold_cutoff),
            arrowPosition="right",
            doShadowText=blockarrow_shadowtext,
            blockWidthPercent=5*blockarrow_cex[1], 
            arrowLabel=label_up, 
            col=blockarrow_colors[["up"]], 
            labelCex=blockarrow_label_cex * blockarrow_cex[1],
            labelFont=blockarrow_font, 
            arrowLabelColor=blockarrow_label_colors[["up"]],
            arrowLabelBorder=alpha2col(blockarrow_label_colors[["up"]], 0.3));
         blockArrowMargin(axisPosition="topAxis", 
            xright=-log2(fold_cutoff),
            arrowPosition="left",
            doShadowText=blockarrow_shadowtext,
            blockWidthPercent=5*blockarrow_cex[1], 
            arrowLabel=label_dn, 
            col=blockarrow_colors[["down"]], 
            labelCex=blockarrow_label_cex * blockarrow_cex[1],
            labelFont=blockarrow_font, 
            arrowLabelColor=blockarrow_label_colors[["down"]],
            arrowLabelBorder=jamba::alpha2col(blockarrow_label_colors[["down"]], 0.3));
      }

      #parList[["postBlockArrows"]] <- par(no.readonly=TRUE);
      
      ## Display the significance cutoff used
      #subTitle <- paste("Significance cutoff <= ", pvalueCutoff, ", and fold change cutoff > ", round(digits=2, 2^fcCutoff));

      if (do_cutoff_caption) {
         captions <- character(0);
         if (length(sig_cutoff) > 0) {
            sig_cutoff_label <- signif(digits=2, sig_cutoff);
            sig_comp <- "<=";
            if (sig_cutoff_label > sig_cutoff) {
               sig_comp <- "<";
            }
            captions <- c(captions,
               paste(
                  "significance",
                  sig_comp,
                  sig_cutoff_label));
         }
         if (length(fold_cutoff) > 0) {
            fold_cutoff_label <- signif(digits=2, fold_cutoff);
            fold_comp <- ">=";
            if (fold_cutoff_label < fold_cutoff) {
               fold_comp <- ">";
            }
            captions <- c(captions,
               paste(
                  "fold",
                  fold_comp,
                  fold_cutoff_label));
         }
         if (length(expr_cutoff) > 0) {
            expr_cutoff_label <- signif(digits=2, expr_cutoff);
            expr_comp <- ">=";
            if (expr_cutoff_label > expr_cutoff) {
               expr_comp <- ">";
            }
            captions <- c(captions,
               paste(
                  "signal",
                  expr_comp,
                  expr_cutoff_label));
         }
         if (length(captions) == 0) {
            caption <- "No statistical thresholds were applied.";
         } else {
            caption <- paste(captions,
               collapse=", ");
         }
         title(sub=caption,
            adj=0.99,
            cex.sub=caption_cex,
            line=parMar[1]-1.5);
         
         ## Display the total points
         total_sub <- paste0("Total points: ",
            format(big.mark=",", nrow(x)));
         title(sub=total_sub,
            adj=0.01,
            cex.sub=caption_cex,
            line=parMar[1] - 1.5);
      }
      
      ## Overall Title
      if (!tophist) {
         do_overall_title();
      }
   }
   
   
   if (1 == 2 && tophist) {
      par("mar"=origPar$mar);
      par("plt"=origPar$plt);
      par("usr"=origPar$usr);
   }
   return(invisible(
      list(
         x=x_values,
         y=y_values,
         point_type=point_type,
         point_colors=point_colors,
         border_colors=border_colors,
         hi_points=hi_points,
         parList=parList)));
}

#' Draw block arrows in plot margins
#' 
#' Draw block arrows in plot margins
#' 
#' This function draws block arrows in plot margins,
#' intended to be used to describe plot axis ranges
#' with an optional label. The driving example is to
#' describe the number of statistical hits shown
#' on a volcano plot.
#' 
#' Run `blockArrowMargin(doExample=TRUE)` to see a visual
#' example.
#' 
#' @examples
#' blockArrowMargin(doExample=TRUE)
#' 
#' @export
blockArrowMargin <- function
(axisPosition="rightAxis",
 arrowPosition="top",
 arrowLabel="",
 arrowDirection="updown",
 xleft=NULL,
 xright=NULL,
 ybottom=NULL,
 ytop=NULL,
 col="#660000FF",
 labelFont=1,
 labelCex=1,
 border="#000000FF",
 xpd=TRUE,
 parUsr=par("usr"),
 arrowWidthPercent=6,
 blockWidthPercent=5,
 arrowLengthPercent=blockWidthPercent*0.6,
 bufferPercent=0.5,
 blankFirst=FALSE,
 arrowLabelColor=NULL,#"#FFFFFFFF",
 arrowLabelBorder="#000000FF",
 doExample=FALSE,
 doBlockGradient=TRUE,
 doShadowText=TRUE,
 gradientDarkFactor=1.5,
 gradientSFactor=1.5,
 verbose=FALSE,
 ...)
{
   ## Purpose is to draw a block arrow, like ones you see in PowerPoint, but using
   ## rect() syntax as if drawing a rectangle.
   ##
   ## Currently the function draws block arrows outside the plot, as if to label
   ## an axis.
   ##
   ## Trim some of the width away to give room for the arrow to be drawn.
   ##
   ## Some examples:
   if (doShadowText) {
      text_fun <- jamba::shadowText;
   } else {
      text_fun <- text;
   }
   
   if (doExample) {
      jamba::nullPlot();
      blockArrowMargin(axisPosition="rightAxis",
         ybottom=1.52, 
         arrowPosition="top", 
         col="#BB0000FF", 
         arrowLabel="Up-regulated",
         arrowWidthPercent=arrowWidthPercent, 
         arrowLengthPercent=arrowLengthPercent,
         blockWidthPercent=blockWidthPercent, 
         labelCex=labelCex, 
         ...);
      blockArrowMargin(axisPosition="rightAxis", 
         ytop=1.48, 
         arrowPosition="bottom", 
         col="#0000BBFF", 
         arrowLabel="Down-regulated",
         arrowWidthPercent=arrowWidthPercent, 
         arrowLengthPercent=arrowLengthPercent,
         blockWidthPercent=blockWidthPercent, 
         labelCex=labelCex, 
         ...);
      blockArrowMargin(axisPosition="topAxis", 
         xright=1.48, 
         arrowPosition="left", 
         col="#0000BBFF", 
         arrowLabel="Down-regulated",
         arrowWidthPercent=arrowWidthPercent,
         arrowLengthPercent=arrowLengthPercent,
         blockWidthPercent=blockWidthPercent, 
         labelCex=labelCex, 
         ...);
      blockArrowMargin(axisPosition="topAxis", 
         xleft=1.52, 
         arrowPosition="right", 
         col="#BB0000FF", 
         arrowLabel="Up-regulated",
         arrowWidthPercent=arrowWidthPercent, 
         arrowLengthPercent=arrowLengthPercent,
         blockWidthPercent=blockWidthPercent, 
         labelCex=labelCex, 
         ...);
      return(NULL);
   }
   if (length(axisPosition) > 1) {
      retVals <- lapply(axisPosition, function(iAxisPosition) {
         blockArrowMargin(axisPosition=iAxisPosition, 
            arrowPosition=arrowPosition, 
            arrowLabel=arrowLabel,
            arrowDirection=arrowDirection, 
            xleft=xleft, 
            xright=xright, 
            ybottom=ybottom, 
            ytop=ytop,
            col=col, 
            labelFont=1, 
            labelCex=1, 
            border=border, 
            xpd=xpd, 
            parUsr=parUsr,
            arrowWidthPercent=arrowWidthPercent, 
            arrowLengthPercent=arrowLengthPercent,
            blockWidthPercent=blockWidthPercent, 
            bufferPercent=bufferPercent,
            blankFirst=blankFirst,
            arrowLabelColor=arrowLabelColor, 
            arrowLabelBorder=arrowLabelBorder,
            doExample=doExample, 
            doBlockGradient=doBlockGradient, 
            gradientDarkFactor=gradientDarkFactor, 
            gradientSFactor=gradientSFactor, 
            ...);
      });
      invisible(retVals);
   }
   
   ## figure aspect, greater than 1 is wider than tall
   parFin <- par("fin");
   fig_aspect <- par("fin")[1] / par("fin")[2];
   if (jamba::igrepHas("left|right", axisPosition)) {
      if (fig_aspect < 1) {
         blockWidthPercent <- blockWidthPercent * fig_aspect;
         arrowLengthPercent <- arrowLengthPercent * 1;
      } else {
         blockWidthPercent <- blockWidthPercent * 1;
         arrowLengthPercent <- arrowLengthPercent * fig_aspect;
      }
   }
   if (jamba::igrepHas("top|bottom", axisPosition)) {
      if (fig_aspect > 1) {
         blockWidthPercent <- blockWidthPercent * fig_aspect;
         arrowLengthPercent <- arrowLengthPercent / fig_aspect;
      } else {
         blockWidthPercent <- blockWidthPercent * 1;
         arrowLengthPercent <- arrowLengthPercent * fig_aspect;
      }
   }

   
   if (doBlockGradient) {
      col <- jamba::alpha2col(
         alpha=jamba::col2alpha(col),
         jamba::fixYellow(col));
      if (length(col) == 1) {
         col2 <- jamba::fixYellow(
            jamba::makeColorDarker(
               col,
               darkFactor=gradientDarkFactor, 
               sFactor=gradientSFactor));
         col2 <- jamba::alpha2col(col2,
            alpha=jamba::col2alpha(col));
         if (jamba::igrepHas("topbottom|bottomtop|leftright|rightleft", arrowPosition)) {
            col <- c(col2, col, col, col2);
         } else {
            col <- c(col, col2);
         }
         colGradient <- colorRampPalette(col,
            alpha=TRUE)(35);
      } else {
         colGradient <- colorRampPalette(col,
            alpha=TRUE)(35);
         col2 <- tail(colGradient, 1);
         col <- head(colGradient, 1);
      }
      col1 <- col[1];
      if (length(arrowLabelColor) == 0) {
         arrowLabelColor <- jamba::setTextContrastColor(head(col2, 1));
      }
   }
   if (length(arrowLabelColor) == 0) {
      arrowLabelColor <- jamba::setTextContrastColor(head(col, 1));
   }

   arrowSets <- list("empty"=list(x=NULL, y=NULL));
   plotWidth <- diff(parUsr[1:2]) / diff(par("plt")[1:2]);
   plotHeight <- diff(parUsr[3:4]) / diff(par("plt")[3:4]);
   if (jamba::igrepHas("rightAxis|leftAxis", axisPosition)) {
      if (is.null(ybottom)) {
         ybottom <- parUsr[3];
      }
      if (is.null(ytop)) {
         ytop <- parUsr[4];
      }
      if (jamba::igrepHas("rightAxis", axisPosition)) {
         if (is.null(xleft)) {
            xleft <- parUsr[2] + 
               plotWidth * bufferPercent / 100;
         }
         if (is.null(xright)) {
            xright <- parUsr[2] + 
               plotWidth * blockWidthPercent/100 +
               plotWidth * bufferPercent / 100;
         }
      } else {
         if (is.null(xright)) {
            xright <- parUsr[1] -
               plotWidth * bufferPercent / 100;
         }
         if (is.null(xleft)) {
            xleft <- parUsr[1] -
               plotWidth * blockWidthPercent/100 -
               plotWidth * bufferPercent / 100;
         }
      }
      arrowDirection <- "updown";
      srtLabel <- 90;
      gradientXY <- "y";
   } else if (jamba::igrepHas("topAxis|bottomAxis", axisPosition)) {
      if (jamba::igrepHas("topAxis", axisPosition)) {
         if (is.null(ybottom)) {
            ybottom <- parUsr[4] +
               plotHeight * bufferPercent / 100;
         }
         if (is.null(ytop)) {
            ytop <- parUsr[4] +
               plotHeight * blockWidthPercent/100 +
               plotHeight * bufferPercent / 100;
         }
      } else {
         if (is.null(ytop)) {
            ytop <- parUsr[3] -
               plotHeight * bufferPercent / 100;
         }
         if (is.null(ybottom)) {
            ybottom <- parUsr[3] -
               plotHeight * blockWidthPercent/100 -
               plotHeight * bufferPercent / 100;
         }
      }
      if (is.null(xleft)) {
         xleft <- parUsr[1];
      }
      if (is.null(xright)) {
         xright <- parUsr[2];
      }
      arrowDirection <- "leftright";
      srtLabel <- 0;
      gradientXY <- "x";
   }
   if (blankFirst) {
      #printDebug("blanking the area first");
      polygon(x=c(xleft, xright, xright, xleft),
         y=c(ybottom, ybottom, ytop, ytop), 
         col="white", 
         border="white", 
         xpd=TRUE);
   }
   if (length(jamba::igrep("up|down", arrowDirection)) > 0) {
      ## Trim away the y-coordinates
      arrowWidth <- abs(xright - xleft);
      arrowWidthDiff <- arrowWidth * (arrowWidthPercent / 50);
      if (xright > xleft) {
         xright1 <- xright - arrowWidthDiff;
         xleft1 <- xleft + arrowWidthDiff;
      } else {
         xright1 <- xright + arrowWidthDiff;
         xleft1 <- xleft - arrowWidthDiff;
      }
      arrowLength <- abs(ytop - ybottom);
      arrowLengthDiff <- plotHeight * arrowLengthPercent / 100;
      if (arrowLengthDiff >= arrowLength) {
         arrowLengthDiff <- arrowLength / 3;
      }
      if (ytop > ybottom) {
         if (length(jamba::igrep("top", arrowPosition)) > 0) {
            ytop1 <- ytop - arrowLengthDiff;
            yArrowPoints <- c(ytop1, ytop1, ytop, ytop1, ytop1);
            xArrowPoints <- c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop1, ybottom));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ytop, ytop);
            xArrowPoints <- c(xleft1, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
         if (length(jamba::igrep("bottom", arrowPosition)) > 0) {
            ybottom1 <- ybottom + arrowLengthDiff;
            yArrowPoints <- c(ybottom1, ybottom1, ybottom, ybottom1, ybottom1);
            xArrowPoints <- rev(c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1));
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom1));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ybottom, ybottom);
            xArrowPoints <- c(xright1, xleft1);
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
      } else {
         if (length(jamba::igrep("bottom", arrowPosition)) > 0) {
            ytop1 <- ytop + arrowLengthDiff;
            yArrowPoints <- c(ytop1, ytop1, ytop, ytop1, ytop1);
            xArrowPoints <- rev(c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1));
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop1, ybottom));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ytop, ytop);
            xArrowPoints <- c(xright1, xleft1);
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
         if (length(jamba::igrep("top", arrowPosition)) > 0) {
            ybottom1 <- ybottom - arrowLengthDiff;
            yArrowPoints <- c(ybottom1, ybottom1, ybottom, ybottom1, ybottom1);
            xArrowPoints <- c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom1));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ybottom, ybottom);
            xArrowPoints <- c(xleft1, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
      }
   } else {
      ## Arrow position is left-to-right
      ## Trim away the x-coordinates
      arrowWidth <- abs(ytop - ybottom);
      arrowWidthDiff <- arrowWidth * (arrowWidthPercent / 50);
      if (ytop > ybottom) {
         ytop1 <- ytop - arrowWidthDiff;
         ybottom1 <- ybottom + arrowWidthDiff;
      } else {
         ytop1 <- ytop + arrowWidthDiff;
         ybottom1 <- ybottom - arrowWidthDiff;
      }
      arrowLength <- abs(xright - xleft);
      arrowLengthDiff <- plotWidth * arrowLengthPercent / 100;
      if (arrowLengthDiff >= arrowLength) {
         arrowLengthDiff <- arrowLength / 3;
      }

      if (xright > xleft) {
         if (length(jamba::igrep("right", arrowPosition)) > 0) {
            xright1 <- xright - arrowLengthDiff;
            xArrowPoints <- c(xright1, xright1, xright, xright1, xright1);
            yArrowPoints <- c(ybottom1, ybottom, mean(c(ybottom, ytop)), ytop, ytop1);
            arrowSets <- c(arrowSets, list("right"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright1));
         } else {
            xright1 <- xright;
            xArrowPoints <- c(xright, xright);
            yArrowPoints <- c(ybottom1, ytop1);
            arrowSets <- c(arrowSets, list("right"=list(x=xArrowPoints, y=yArrowPoints)));
            #print(arrowSets);
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
         if (length(jamba::igrep("left", arrowPosition)) > 0) {
            xleft1 <- xleft + arrowLengthDiff;
            xArrowPoints <- c(xleft1, xleft1, xleft, xleft1, xleft1);
            yArrowPoints <- rev(c(ybottom1, ybottom, mean(c(ybottom, ytop)), ytop, ytop1));
            arrowSets <- c(arrowSets, list("left"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft1, xright1));
         } else {
            xleft1 <- xleft;
            xArrowPoints <- c(xleft, xleft);
            yArrowPoints <- rev(c(ybottom1, ytop1));
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft1, xright1));
         }
      } else {
      }
      
   }
   yLabel <- mean(c(ytop, ybottom));
   xLabel <- mean(c(xleft, xright));
   
   allX <- unlist(lapply(jamba::unvigrep("^empty$", names(arrowSets)), function(i){
      j <- arrowSets[[i]];
      j["x"];
   }))
   allY <- unlist(lapply(jamba::unvigrep("^empty$", names(arrowSets)), function(i){
      j <- arrowSets[[i]];
      j["y"];
   }));
   if (doBlockGradient) {
      arrowSides <- sapply(jamba::unvigrep("^empty$", names(arrowSets)), function(i){
         length(arrowSets[[i]][["x"]]) > 2;
      });
      arrowSides <- paste(names(arrowSides)[arrowSides], collapse="");
      if (arrowSides %in% c("bottom", "left")) {
         colGradient <- rev(colGradient);
         col21 <- col1;
         col1 <- col2;
         col2 <- col1;
      }
      arrowBox <- lapply(jamba::nameVector(jamba::unvigrep("^empty$", names(arrowSets))), function(i){
         xi <- arrowSets[[i]][["x"]];
         yi <- arrowSets[[i]][["y"]];
         xI <- unique(c(xi[1], xi[length(xi)]));
         yI <- unique(c(yi[1], yi[length(yi)]));
         list(x=xI, y=yI);
      });
      arrowBoxX <- unique(sort(c(arrowBox[[1]][["x"]], arrowBox[[2]][["x"]])));
      arrowBoxY <- unique(sort(c(arrowBox[[1]][["y"]], arrowBox[[2]][["y"]])));
      
      xpdPar <- par("xpd");
      par("xpd"=TRUE);
      polygon(x=allX, 
         y=allY, 
         col=tail(colGradient,1), 
         border=border, 
         xpd=TRUE);
      gr1 <- gradient_rect(col=colGradient, 
         gradient=gradientXY, 
         xleft=arrowBoxX[1], 
         xright=arrowBoxX[2],
         ybottom=arrowBoxY[1], 
         ytop=arrowBoxY[2], 
         border="#00000000");
      par("xpd"=xpdPar);
      arrowsDrawn1 <- lapply(jamba::nameVector(jamba::vgrep("top|right", names(arrowSets))), function(i){
         as1 <- arrowSets[[i]];
         polygon(x=as1[["x"]], 
            y=as1[["y"]], 
            col=col2, 
            border=NA, 
            xpd=TRUE);
         as1;
      });
      arrowsDrawn2 <- lapply(jamba::nameVector(jamba::vgrep("bottom|left", names(arrowSets))), function(i){
         as1 <- arrowSets[[i]];
         polygon(x=as1[["x"]], 
            y=as1[["y"]], 
            col=col1, 
            border=col1, 
            xpd=TRUE);
         as1;
      });
      polygon(x=allX, 
         y=allY, 
         col=NA, 
         border=border, 
         xpd=TRUE);
   } else {
      polygon(x=allX, 
         y=allY, 
         col=col, 
         border=border, 
         xpd=TRUE);
      arrowsDrawn1 <- NULL;
      arrowsDrawn2 <- NULL;
   }
   
   ## Optionally label the block arrows
   if (!is.null(arrowLabel) && !arrowLabel %in% c(NA, "")) {
      if (verbose) {
         jamba::printDebug("xLabel: ", 
            round(digits=2, xLabel), 
            ",  yLabel: ", 
            round(digits=2, yLabel));
      }
      text_fun(label=arrowLabel,
         x=xLabel,
         y=yLabel,
         srt=srtLabel,
         xpd=TRUE,
         col=arrowLabelColor,
         font=labelFont,
         cex=labelCex,
         #bg=arrowLabelBorder,
         adj=c(0.5,0.5), 
         ...);
   }
   retVals <- list(arrowSets=arrowSets, 
      arrowsDrawn1=arrowsDrawn1, 
      arrowsDrawn2=arrowsDrawn2,
      allX=allX, 
      allY=allY);
   invisible(retVals);
}

#' Rectangle with color gradient fill
#' 
#' Rectangle with color gradient fill
#' 
#' This function was inspired by the `plotrix::gradient.rect()`
#' function in the plotrix R package. The function is
#' simplified here, and requires a vector of colors in `col`.
#' 
#' @param xleft,ybottom,xright,ytop `numeric` vectors indicating
#'    the position of sides of a rectangle, passed to
#'    `graphics::rect()`. Multiple rectangles may be defined.
#' @param col `character` vector of colors used to fill the rectangles.
#' @param gradient `character` string indicating the direction of
#'    color gradient, with two allowed values: `"x"` and `"y"`.
#' @param border `character` value indicating the color of border
#'    around the rectangle.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' jamba::nullPlot(xlim=c(0,5), ylim=c(0,5), xaxt="s", yaxt="s");
#' gradient_rect(xleft=1, 
#'    ybottom=1, 
#'    xright=2.5, 
#'    ytop=2.5, 
#'    col=jamba::getColorRamp("Reds", n=15))
#' gradient_rect(xleft=2.5, 
#'    ybottom=2.5, 
#'    xright=4, 
#'    ytop=4, 
#'    gradient="y",
#'    col=jamba::getColorRamp("Reds", n=15))
#' 
#' @export
gradient_rect <- function
(xleft,
 ybottom,
 xright,
 ytop,
 col,
 gradient="x",
 border=par("fg"),
 ...)
{
   nslices <- length(col)

   nrect <- max(unlist(lapply(list(xleft, ybottom, xright, ytop),
      length)));
   oldxpd <- par(xpd = NA)
   if (nrect > 1) {
      if (length(xleft) < nrect)
         xleft <- rep(xleft, length.out=nrect)
      if (length(ybottom) < nrect)
         ybottom <- rep(ybottom, length.out=nrect)
      if (length(xright) < nrect)
         xright <- rep(xright, length.out=nrect)
      if (length(ytop) < nrect)
         ytop <- rep(ytop, length.out=nrect)
      for (i in 1:nrect) {
         gradient_rect(xleft[i],
            ybottom[i],
            xright[i],
            ytop[i],
            col=col,
            nslices=nslices,
            gradient=gradient,
            border=border,
            ...)
      }
   } else {
      if (gradient == "x") {
         xinc <- (xright - xleft)/nslices;
         xlefts <- seq(xleft,
            xright - xinc,
            length=nslices);
         xrights <- xlefts + xinc;
         rect(xlefts,
            ybottom,
            xrights,
            ytop,
            col=col,
            lty=0);
         rect(xlefts[1],
            ybottom,
            xrights[nslices],
            ytop,
            border=border);
      } else {
         yinc <- (ytop - ybottom)/nslices;
         ybottoms <- seq(ybottom,
            ytop - yinc,
            length=nslices);
         ytops <- ybottoms + yinc;
         rect(xleft,
            ybottoms,
            xright,
            ytops,
            col=col,
            lty=0);
         rect(xleft,
            ybottoms[1],
            xright,
            ytops[nslices],
            border=border);
      }
   }
   par(oldxpd);
   invisible(col);
}


#' Log-scaled axis including transformed P-values
#' 
#' @export
logAxis <- function
(side,
 at=NULL,
 base=2,
 values=FALSE,
 useFcValues=TRUE,
 las=2,
 makeNegative=FALSE,
 doSignedSignificance=FALSE,
 flipSignedSignificanceSign=FALSE,
 bigMark=",",
 prettyN=5,
 digits=2,
 cex.axis=1,
 font.axis=1,
 ...)
{
   ## Purpose is to draw an axis using log scale
   ## Logic borrowed from 'log10' package, but extended
   ## to allow 2-based (or N-based) log transforms
   if (tolower(side) %in% c("x", "bottom")) {
      side <- 1;
   }
   if (tolower(side) %in% c("y", "left")) {
      side <- 2;
   }
   if (tolower(side) %in% c("above", "top")) {
      side <- 3;
   }
   if (tolower(side) %in% c("right")) {
      side <- 4;
   }
   if (is.null(at)) {
      if (side %in% c(1,3)) {
         #at1 <- rmNA(log(pretty(base^par("usr")[1:2], n=prettyN*10), base=base));
         at <- pretty(c(par("usr")[1:2]), n=prettyN);
      } else {
         #at1 <- rmNA(log(pretty(base^par("usr")[3:4], n=prettyN*10), base=base));
         at <- pretty(c(par("usr")[3:4]), n=prettyN);
      }
      #printDebug(c("at: ", paste(at, collapse=", ")));
   }
   sa1 <- lapply(at, function(i){
      b <- as.numeric(base);
      j <- as.numeric(i);
      #doSignedSignificance <<- doSignedSignificance;
      if (doSignedSignificance %in% c(TRUE, "TRUE")) {
         ## For this operation, force not displaying the value itself
         values <- FALSE;
         ## Grab the sign from the exponent
         jSign <- sign(j);
         ## Make exponents always negative
         j1 <- -(abs(j));
         ## Optionally flip the sign
         if (flipSignedSignificanceSign) {
            b <- jSign * abs(b);
         }
      } else if (makeNegative) {
         j1 <- -j;
      } else {
         j1 <- j;
      }
      if (values) {
         xLabels <- b^j1;
         ## For 2^(-2), instead of reporting 0.5, report -2
         if (useFcValues) {
            xLabels[xLabels < 1] <- (-1 / xLabels[xLabels < 1]);
            xLabels <- signif(xLabels, digits=digits);
         }
         if (!is.null(bigMark) && !bigMark %in% c("")) {
            if (abs(xLabels) >= 1) {
               xLabels <- signif(xLabels, digits=digits);
            }
            xLabels <- format(xLabels, big.mark=bigMark, trim=TRUE);
         }
         axis(side=side,
            at=j,
            labels=xLabels,
            las=2,
            cex.axis=cex.axis,
            font.axis=font.axis,
            ...);
      } else {
         if (i == 0) {
            axis(side=side,
               at=j,
               labels=1,
               las=2,
               cex.axis=cex.axis*1,
               font.axis=font.axis,
               ...);
            xLabels <- 1;
         } else {
            xLabels <- b^j1;
            b <- as.character(b);
            j1 <- as.character(j1);
            if (font.axis == 1) {
               if (grepl("e", format(xLabels))) {
                  # if format() would use exponent, we display exponent
                  axis(side=side,
                     at=j,
                     labels=substitute(b^j1),
                     las=2,
                     cex.axis=cex.axis*1,
                     font.axis=font.axis,
                     ...);
               } else {
                  # if format() would not use exponent, we do not display exponent
                  axis(side=side,
                     at=j,
                     labels=format(xLabels),
                     las=2,
                     cex.axis=cex.axis*1,
                     font.axis=font.axis,
                     ...);
               }
            } else {
               str <- paste0('axis(side, at=j, labels=expression(bold("', b, '"^"', j1, '")), las=2, cex.axis=cex.axis*1, font.axis=2)')
               eval(parse(text=str));
            }
            xLabels <- substitute(b^j1);
         }
      }
      list(j=j, j1=j1, b=b, xLabels=xLabels);
   });
   invisible(list(sa1));
}

