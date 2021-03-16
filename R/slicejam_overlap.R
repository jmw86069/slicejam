#
# slicejam_overlap.R
#
# functions that analyze peak overlaps


#' Peak overlap calculations
#' 
#' Peak overlap calculations
#' 
#' This function calculates all pairwise overlap counts for
#' a list of peak sets.
#' 
#' Note that the number of overlapping peaks is not
#' always symmetric, for example 100 peaks in `SetA`
#' might overlap 75 peaks in `SetB`.
#' 
#' The function can optionally filter peaks using `min_score`,
#' which is applied to each set of peaks in `grl`. Note that
#' `min_score` can be a vector with length `length(grl)`
#' so that each peak set can be filtered using its own
#' distinct score.
#' 
#' The argument `metric` defines the type of overlap to
#' calculate. Each method can have the suffix `_bases`
#' which means the calculation is based upon the number
#' of base positions, and not the number of peak regions.
#' By default all calculations are performed and returned,
#' however a subset can be provided which may be slightly
#' faster to calculate.
#' 
#' A brief description of terms used in the calculations:
#' 
#' * `num_A`: number of peaks in A
#' * `num_B`: number of peaks in B
#' * `num_A_over_B`: number of peaks in A that overlap any peak in B
#' * `bases_A`: number of unique genome bases covered by peaks in A
#' * `bases_B`: number of unique genome bases covered by peaks in B
#' * `bases_A_over_B`: number of bases where peaks in A overlap peaks in B
#' 
#' A brief description of each metric:
#' 
#' * `count`: num_A_overB
#' * `percent`: num_A_overB / num_A
#' * `maxpercent`: num_A_overB / minimum( num_A, num_B )
#' * `jaccard`: num_A_overB / ( num_A + num_B - num_A_over_B )
#' * `count_bases`: bases_A_overB
#' * `percent_bases`: bases_A_overB / bases_A
#' * `maxpercent_bases`: bases_A_overB / minimum( bases_A, bases_B )
#' * `jaccard_bases`: bases_A_overB / ( bases_A + bases_B - bases_A_over_B )
#' 
#' The argument `count_method` is relevant only to peak-based
#' overlaps, not for the metric `_bases`. It defines whether
#' to return: the number of peaks in A that overlap peaks
#' in B, regardless whether any one peak is overlapped
#' multiple times; or return the total number of overlap
#' events even including peaks with multiple overlaps.
#' In practice the different is usually quite small, and
#' is relevant usually only when the peak sizes between
#' A and B are quite different, or include large regions
#' of maybe 1,000 base width or larger. In these cases,
#' it may be more useful and informative to use the `_bases`
#' metric calculations.
#' 
#' Note that this method has not been tested extensively
#' for strand-specific overlaps, using `ignore.strand=FALSE`,
#' although the calculations metrics should be applicable.
#' One might use strand-specific overlaps only with
#' RNA-seq protocols that are strand-specific, for example
#' to test the number of overlaps with annotated stranded
#' transcript features.
#' 
#' For several of the calculations, the possibility of zero
#' denominator is prevented by forcing minimum value `1`,
#' which would render the result `0` in that case. This
#' result is considered more consistent with saying there is
#' `"nothing"` present for the ratio.
#' 
#' @family slicejam overlap
#' 
#' @return `list` of `matrix` objects that contain overlap
#'    metrics. There will be one entry per input `metric`.
#' 
#' @param grl either `GRangesList`, or `list` of `GRanges`.
#' @param metric `character` string indicating the overlap
#'    metric to calculate. See this function description
#'    for details.
#' @param count_method `string` indicating how to define
#'    overlap counts:
#'    * `"peak"` returns the number of peaks
#'    in `SetA` that overlap one or more peaks in `SetB`;
#'    * `"overlap"` returns the total number of overlap
#'    events between `SetA` and `SetB`.
#' @param min_score `numeric` indicating the minimum score to use
#'    for each set of peaks provided in `grl`. This option requires
#'    the score to be stored in a column `"score"` for each set of peaks.
#'    For any set of peaks that does not contain `"score"` the filter
#'    will not be applied. Note that `min_score` can be a vector with
#'    length `length(grl)` so that each peak set can be filtered using
#'    its own relevant score threshold. When `min_score` contains `NA`
#'    the threshold is not applied. Use `verbose=TRUE` to confirm
#'    this behavior.
#' @param ignore.strand `logical` passed to `GenomicRanges::countOverlaps()`,
#'    `GenomicRanges::reduce()` and `GenomicRanges::intersect()` which
#'    indicates whether strandedness is used for overlap calculations.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... other arguments are ignored.
#' 
#' @examples
#' # Typical workflow starting with GRanges objects in a list
#' # Note: without test data the steps below are not run
#' if (1 == 2) {
#' 
#' # peak overlap calcs
#' groverlaplist <- peakoverlap_calcs(grl)
#' 
#' # you can review one heatmap with other calcs as labels
#' # the easiest is to order the overlap list
#' overlap_order <- c("percent_bases", "percent", "count", "count_bases");
#' overlap_suffix <- c("% bp", "% pk", " pk", " bp")
#' peakoverlap_heatmap(groverlaplist[overlap_order],
#'    suffix=overlap_suffix)
#' 
#' # much fancier method is to put every heatmap into one figure (below)
#' #
#' # generate a list of heatmaps
#' hms <- lapply(jamba::nameVectorN(groverlaplist), function(n){
#'    m <- groverlaplist[[n]];
#'    hm <- peakoverlap_heatmap(m,
#'       column_title=n,
#'          name=n)
#' })
#' 
#' # convert them all to grid objects which makes them "real"
#' hms_grids <- lapply(hms, function(hm){
#'    grid::grid.grabExpr(ComplexHeatmap::draw(hm))
#' })
#' 
#' # save them in a wide array of heatmaps in PDF format
#' outdir <- ".";
#' cairo_pdf(file.path(outdir,
#'    "peakoverlap_heatmaps.pdf"),
#'    onefile=TRUE, width=32, height=14)
#' cowplot::plot_grid(ncol=4,
#'    plotlist=hms_grids[c(1,2,3,4,5,6,7,8)]
#' )
#' dev.off()
#' 
#' }
#' 
#' @export
peakoverlap_calcs <- function
(grl,
 metric=c("count",
    "percent",
    "maxpercent",
    "jaccard",
    "count_bases",
    "percent_bases",
    "maxpercent_bases",
    "jaccard_bases"),
 count_method=c("peak", "overlap"),
 min_score=NULL,
 ignore.strand=TRUE,
 verbose=FALSE,
 ...)
{
   # validate input arguments
   metric <- match.arg(metric,
      several.ok=TRUE);
   count_method <- match.arg(count_method);
   
   # grnames
   if (length(names(grl)) == 0) {
      names(grl) <- seq_along(grl);
   }
   grlnames <- names(grl);
   
   # minimum score per peak set
   if (length(min_score) > 0) {
      min_score <- rep(min_score, length.out=length(grl));
      grl <- lapply(seq_along(grl), function(grnum){
         gr <- grl[[grnum]];
         if (!is.na(min_score[[grnum]]) && "score" %in% colnames(values(gr))) {
            gr <- subset(gr, score >= min_score[[grnum]])
            if (verbose) {
               jamba::printDebug("peakoverlap_calcs(): ",
                  "min_score ", min_score[[grnum]], " applied to ",
                  grlnames[grnum]);
            }
         } else {
            if (verbose) {
               jamba::printDebug("peakoverlap_calcs(): ",
                  "min_score ", min_score[[grnum]], " NOT applied to ",
                  grlnames[grnum],
                  fgText=c("darkorange", "red2"));
            }
         }
         gr
      });
      names(grl) <- grlnames;
   }
   
   # small wrapper method
   get_overlap_counts <- function(a, b, count_method, ignore.strand=TRUE) {
      ij_overlaps <- GenomicRanges::countOverlaps(a,
         b,
         ignore.strand=ignore.strand);
      if ("peak" %in% count_method) {
         overlap <- sum(ij_overlaps > 0);
         nonoverlap <- sum(ij_overlaps == 0);
      } else if ("overlap" %in% count_method) {
         overlap <- sum(ij_overlaps);
         nonoverlap <- sum(ij_overlaps == 0);
      }
      return(list(
         overlap=overlap,
         nonoverlap=nonoverlap));
   }
   
   # overlap counts for all pairwise combinations
   grpairs <- as.list(data.frame(combn(grlnames, 2)));
   if (verbose) {
      jamba::printDebug("peakoverlap_calcs(): ",
         "Calculating ",
         length(grpairs),
         " pairwise overlaps.");
   }
   groverlaps <- jamba::rbindList(lapply(grpairs, function(ij){
      i <- ij[1];
      j <- ij[2];
      ij_overlaps <- get_overlap_counts(grl[[i]], grl[[j]], count_method);
      ijdf <- data.frame(File1=i,
         File2=j,
         overlap=ij_overlaps$overlap,
         nonoverlap=ij_overlaps$nonoverlap)
      ji_overlaps <- get_overlap_counts(grl[[j]], grl[[i]], count_method);
      jidf <- data.frame(File1=j,
         File2=i,
         overlap=ji_overlaps$overlap,
         nonoverlap=ji_overlaps$nonoverlap)
      if (any(grepl("_bases", metric))) {
         # reduce to get base-level coverage
         ij_union <- GenomicRanges::reduce(c(grl[[i]],
            grl[[j]]),
            ignore.strand=TRUE);
         ij_intersect <- GenomicRanges::intersect(grl[[i]],
            grl[[j]],
            ignore.strand=TRUE);
         iwidth <- sum(width(GenomicRanges::reduce(grl[[i]])));
         jwidth <- sum(width(GenomicRanges::reduce(grl[[j]])));
         ijdf$base_overlap <- sum(width(ij_intersect));
         # add columns to ijdf
         ijdf$base_union <- sum(width(ij_union));
         ijdf$base_i <- c(iwidth);
         ijdf$base_j <- c(jwidth);
         # add columns to jidf
         jidf$base_overlap <- sum(width(ij_intersect));
         jidf$base_union <- sum(width(ij_union));
         jidf$base_i <- c(jwidth);
         jidf$base_j <- c(iwidth);
      } else {
      }
      rbind(ijdf, jidf)
   }))
   
   # convert pairwise counts into a numeric matrix
   # percent overlap matrix
   groverlaps2pct <- matrix(numeric(0),
      ncol=length(grlnames),
      nrow=length(grlnames))
   colnames(groverlaps2pct) <- grlnames;
   rownames(groverlaps2pct) <- grlnames;
   if (verbose > 1) {
      jamba::printDebug("peakoverlap_calcs(): ",
         "grlnames:", grlnames);
   }
   # percent bases count matrix
   groverlaps2pctb <- groverlaps2pct;
   # overlap count matrix
   groverlaps2ct <- groverlaps2pct;
   # overlap bases matrix
   groverlaps2ctb <- groverlaps2pct;
   # Jaccard overlap matrix
   groverlaps2j <- groverlaps2pct;
   # Jaccard bases overlap matrix
   groverlaps2jb <- groverlaps2pct;
   # max percent overlap
   groverlaps2max <- groverlaps2pct;
   # max percent bases overlap
   groverlaps2maxb <- groverlaps2pct;
   
   for (k in seq_len(nrow(groverlaps))) {
      i <- groverlaps[k, 1]
      j <- groverlaps[k, 2]
      if (verbose > 1) {
         jamba::printDebug("peakoverlap_calcs(): ",
            "file1:", groverlaps$File1[k],
            ", file2:", groverlaps$File2[k]);
      }
      overlap <- groverlaps[k,"overlap"]
      total <- groverlaps[k,"overlap"] + groverlaps[k,"nonoverlap"]
      # count
      groverlaps2ct[i, j] <- overlap;
      # percent
      if (verbose > 1) {
         jamba::printDebug("peakoverlap_calcs(): ",
            "percent:",
            " overlap:", overlap,
            ", total:", total);
         if (is.na(total)) {
            jamba::printDebug("groverlaps[k,,drop=FALSE]:");
            print(groverlaps[k,,drop=FALSE]);
         }
      }
      groverlaps2pct[i, j] <- (overlap / total);
      # fill the diagonal
      if (is.na(groverlaps2ct[i, i])) {
         groverlaps2ct[i, i] <- length(grl[[i]])
         groverlaps2pct[i, i] <- 1;
      }
      if (is.na(groverlaps2ct[j, j])) {
         groverlaps2ct[j, j] <- length(grl[[j]])
         groverlaps2pct[j, j] <- 1;
      }
      if ("count_bases" %in% metric) {
         if (verbose > 1) {
            jamba::printDebug("peakoverlap_calcs(): ",
               "count_bases");
         }
         # count bases
         groverlaps2ctb[i, j] <- groverlaps[k,"base_overlap"];
         if (is.na(groverlaps2ctb[i, i])) {
            groverlaps2ctb[i, i] <- sum(width(GenomicRanges::reduce(grl[[i]])));
         }
         if (is.na(groverlaps2ctb[j, j])) {
            groverlaps2ctb[j, j] <- sum(width(GenomicRanges::reduce(grl[[j]])));
         }
      }
      if ("percent_bases" %in% metric) {
         if (verbose > 1) {
            jamba::printDebug("peakoverlap_calcs(): ",
               "percent_bases");
         }
         # percent bases
         groverlaps2pctb[i, j] <- groverlaps[k,"base_overlap"] / jamba::noiseFloor(minimum=1, groverlaps[k,"base_i"]);
         if (is.na(groverlaps2pctb[i, i])) {
            groverlaps2pctb[i, i] <- groverlaps[k,"base_i"] / jamba::noiseFloor(minimum=1, groverlaps[k,"base_i"]);
         }
         if (is.na(groverlaps2pctb[j, j])) {
            groverlaps2pctb[j, j] <- groverlaps[k,"base_j"] / jamba::noiseFloor(minimum=1, groverlaps[k,"base_j"]);
         }
      }
      if ("jaccard_bases" %in% metric) {
         if (verbose > 1) {
            jamba::printDebug("peakoverlap_calcs(): ",
               "jaccard_bases");
         }
         # jaccard bases
         groverlaps2jb[i, j] <- groverlaps[k,"base_overlap"] / jamba::noiseFloor(minimum=1, groverlaps[k,"base_union"]);
         groverlaps2jb[i, i] <- 1;
         groverlaps2jb[j, j] <- 1;
      }
      if ("maxpercent_bases" %in% metric) {
         if (verbose > 1) {
            jamba::printDebug("peakoverlap_calcs(): ",
               "maxpercent_bases");
            print(groverlaps[k,,drop=FALSE])
         }
         # maxpercent bases
         base_min <- min(c(groverlaps[k,"base_i"], groverlaps[k,"base_j"]));
         if (verbose > 1) {
            jamba::printDebug("peakoverlap_calcs(): ",
               "base_min:", base_min);
            jamba::printDebug("peakoverlap_calcs(): ",
               "groverlaps2maxb:");
            print(groverlaps2maxb)
         }
         groverlaps2maxb[i, j] <- groverlaps[k,"base_overlap"] / jamba::noiseFloor(minimum=1, base_min);
         groverlaps2maxb[i, i] <- groverlaps[k,"base_i"] / jamba::noiseFloor(minimum=1, groverlaps[k,"base_i"]);
         groverlaps2maxb[j, j] <- groverlaps[k,"base_j"] / jamba::noiseFloor(minimum=1, groverlaps[k,"base_j"]);
      }
   }
   
   # optional
   # adjust diagonal to be pct max peak count
   #diag(groverlaps2m) <- diag(groverlaps2ct) / max(diag(groverlaps2ct))
   
   ## These calculations work best once the full matrix is filled
   #
   # Jaccard overlap coefficient
   # overlap / (sizeA + sizeB - overlap)
   if (verbose > 1) {
      jamba::printDebug("peakoverlap_calcs(): ",
         "jaccard/maxpercent");
   }
   for (i in grlnames) {
      for (j in grlnames) {
         overlap <- groverlaps2ct[i,j];
         union <- groverlaps2ct[i,i] + groverlaps2ct[j,j] - overlap;
         minsize <- min(c(groverlaps2ct[i,i], groverlaps2ct[j,j]));
         # Jaccard overlap
         groverlaps2j[i,j] <- overlap / union;
         # max percent overlap
         groverlaps2max[i,j] <- overlap / minsize;
      }
   }
   
   retlist <- list();
   if ("count" %in% metric) {
      retlist$count <- groverlaps2ct;
   }
   if ("percent" %in% metric) {
      retlist$percent <- groverlaps2pct * 100;
   }
   if ("maxpercent" %in% metric) {
      retlist$maxpercent <- groverlaps2max * 100;
   }
   if ("jaccard" %in% metric) {
      retlist$jaccard <- groverlaps2j;
   }
   if ("count_bases" %in% metric) {
      retlist$count_bases <- groverlaps2ctb;
   }
   if ("percent_bases" %in% metric) {
      retlist$percent_bases <- groverlaps2pctb * 100;
   }
   if ("maxpercent_bases" %in% metric) {
      retlist$maxpercent_bases <- groverlaps2maxb * 100;
   }
   if ("jaccard_bases" %in% metric) {
      retlist$jaccard_bases <- groverlaps2jb;
   }
   return(retlist);
}


#' Heatmap of peak overlap statistics
#' 
#' Heatmap of peak overlap statistics
#' 
#' This function is a method to create a heatmap for the output
#' of `peakoverlap_calcs()`.
#' 
#' See `peakoverlap_calcs()` for the full example workflow.
#' 
#' @family slicejam overlap
#' 
#' @return `ComplexHeatmap::Heatmap` object. Use `ComplexHeatmap::draw()`
#'    or simply `print()` to display the heatmap.
#' 
#' @param m `numeric` matrix
#' @param col any input compatible with `colorjam::getColorRamp()` which
#'    includes: `character` string of a specific color gradient
#'    from `RColorBrewer` or `viridis`; `character` vector of
#'    R colors.
#' @param row_split `NULL`, or vector of `character` or `factor` values
#'    passed to `ComplexHeatmap::Heatmap()` argument `row_split` and
#'    `column_split` to split rows and columns into groups.
#' @param breaks `numeric` vector that defines color breaks used in
#'    the color mapping function `circlize::colorRamp2()`.
#' @param show_breaks `numeric` vector the defines color break positions
#'    to display in the color legend. When `show_breaks=NULL` it uses
#'    `breaks` by default.
#' @param col_hm `NULL` or `function` generated by `circlize::colorRamp2()`
#'    that maps specific numeric ranges to specific colors. When `col_hm`
#'    is supplied, `col` is ignored.
#' @param lens `numeric` value passed to `colorjam::getColorRamp()` which
#'    defines any compression or expansion of the color gradient, where
#'    values higher than 1 make colors more intense, and values less than
#'    -1 make colors less intense across the color range.
#' @param cluster_rows,cluster_columns `logical` indicating whether
#'    to cluster and thus re-order rows and columns.
#' @param m2 `matrix` or `list` of `matrix` objects. When supplied,
#'    each `matrix` is used to add numeric labels to each cell,
#'    relevant only when `show_label=TRUE`.
#' @param prefix,suffix `character` vector of prefix and suffix
#'    values to add to each matrix value displayed when
#'    `show_label=TRUE`. The purpose is to add a label describing the
#'    type of value, for example `suffix="j"` may help reinforce that
#'    a numeric value is a Jaccard overlap coefficient; `suffix="%"`
#'    may indicate a percentage. Note that these labels are applied
#'    to `m` and each matrix in `m2` in order.
#' @param show_label `logical` indicating whether to define
#'    `cell_fun=cell_fun_label()` to display cell labels in the heatmap.
#' @param label_cex `numeric` value passed to `cell_fun_label()` to
#'    adjust the fontsize in for cell labels, as relevant.
#' @param abbrev `logical` passed to `cell_fun_label()` as needed,
#'    which determines whether the numeric value should be abbreviated
#'    for every `1,000` unit of measure. For example 1,500 will be
#'    labeled `1.5k`, and 21,300,000 will be labeled `21.3M`.
#' @param cell_fun `function` or when `cell_fun=NULL` and `show_label=TRUE`
#'    a custom function `cell_fun_label()` will be used that displays the
#'    numeric value in `m` for each cell. If `m2` is also supplied then the
#'    value from each `matrix` in `m2` will also be displayed in each
#'    cell. To force no `cell_fun` to be used, set `cell_fun=FALSE`, or
#'    set `show_label=FALSE` and `cell_fun=NULL`.
#' @param cell_outline `logical` passed to `cell_fun_label()` which is
#'    used only when `show_label=TRUE`. This argument defines whether an
#'    outline is drawn around each heatmap cell.
#' 
#' @export
peakoverlap_heatmap <- function
(m,
   col="Reds",
   row_split=NULL,
   breaks=NULL,
   show_breaks=NULL,
   col_hm=NULL,
   lens=2,
   cluster_rows=FALSE,
   cluster_columns=FALSE,
   m2=NULL,
   prefix="",
   suffix="",
   show_label=TRUE,
   label_cex=1,
   abbrev=TRUE,
   cell_fun=NULL,
   cell_outline=FALSE,
   ...)
{
   # m as a list of matrices
   if (is.list(m)) {
      if (length(m) > 0) {
         m2 <- m[-1];
      }
      m <- m[[1]];
   }
   
   # define color function breaks
   if (length(col_hm) == 0) {
      # define color breaks
      if (length(breaks) < 2) {
         maxm <- max(m, na.rm=TRUE);
         if (maxm <= 1) {
            breaks <- seq(from=0,
               to=1,
               length.out=11)
         } else {
            breaks <- pretty(
               c(0, maxm),
               n=10);
            breaks <- breaks[breaks <= maxm & breaks >= 0];
            breaks1 <- seq(from=0,
               to=max(m, na.rm=TRUE),
               length.out=11)
         }
      }
      if (length(show_breaks) == 0) {
         show_breaks <- breaks;
      }
      col_hm <- circlize::colorRamp2(
         colors=getColorRamp(col,
            lens=lens,
            n=length(breaks)),
         breaks=breaks)
   }
   
   # cell labels
   if (isFALSE(cell_fun)) {
      cell_fun <- NULL;
   } else if (show_label && length(cell_fun) == 0) {
      cell_fun <- cell_fun_label(
         m=c(list(m=m), m2),
         prefix=prefix,
         suffix=suffix,
         cex=label_cex,
         outline=cell_outline,
         col_hm=col_hm,
         abbrev=abbrev,
         ...)
   }
   
   # heatmap
   hm_overlap <- tryCatch({
      ComplexHeatmap::Heatmap(
         m,
         border=TRUE,
         column_split=row_split,
         row_split=row_split,
         cluster_rows=cluster_rows,
         cluster_columns=cluster_columns,
         cell_fun=cell_fun,
         heatmap_legend_param=list(
            at=show_breaks,
            border=TRUE,
            color_bar="discrete"),
         col=col_hm,
         ...)
   }, error=function(e){
      ComplexHeatmap::Heatmap(
         m,
         border=TRUE,
         column_split=row_split,
         row_split=row_split,
         cluster_rows=cluster_rows,
         cluster_columns=cluster_columns,
         cell_fun=cell_fun,
         heatmap_legend_param=list(
            at=show_breaks,
            border=TRUE,
            color_bar="discrete"),
         col=col_hm)
   })
   return(hm_overlap)
}


#' ComplexHeatmap cell function label
#' 
#' ComplexHeatmap cell function label
#' 
#' This function serves as a convenient method to add text
#' labels to each cell in a heatmap produced by
#' `ComplexHeatmap::Heatmap()`, via the argument `cell_fun`.
#' 
#' Note that in its current form, this function requires
#' having the specific color function used for the heatmap,
#' because it uses that color to choose an appropriate
#' contrast color for each cell label. Currently this
#' step requires that the first label matrix is numeric,
#' whose values can be passed to `col_hm` which is the
#' color function used for the heatmap.
#' 
#' This function is slightly unique in that it allows multiple
#' labels, if `m` is supplied as a `list` of `matrix` objects.
#' In fact, some `matrix` objects may contain `character`
#' values with custom labels.
#' 
#' Cell labels are colored based upon the heatmap cell color,
#' which is passed to `jamba::setTextContrastColor()` to determine
#' whether to use light or dark text color for optimum contrast.
#' 
#' Todo: Enable some matrix values that contain `character` data
#' to use `gridtext` for custom markdown formatting. That process
#' requires a slightly different method.
#' 
#' @param m `numeric` matrix or `list` of `matrix` objects. The
#'    first `matrix` object must be `numeric` and compatible
#'    with the color function `col_hm`.
#' @param prefix,suffix `character` vectors that define a prefix and
#'    suffix for each value in `m` for each cell.
#' @param cex `numeric` adjustment for the fontsize used for each label
#' @param col_hm `function` as returned by `circlize::colorRamp2()` which
#'    should be the same function used to create the heatmap
#' @param outline `logical` indicating whether to draw an outline around
#'    each heatmap cell
#' @param abbrev `logical` indicating whether numeric values should
#'    be abbreviated using `jamba::asSize(..., kiloSize=1000)` which
#'    effectively reduces large numbers to `k` for thousands, `M` for
#'    millions (M for Mega), `G` for billions (G for Giga), etc.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' m <- matrix(rnorm(16)*2, ncol=4)
#' colnames(m) <- LETTERS[1:4]
#' rownames(m) <- letters[1:4]
#' col_hm <- circlize::colorRamp2(breaks=(-2:2) * 2,
#'    colors=c("navy", "dodgerblue", "white", "tomato", "red4"))
#' 
#' # the heatmap can be created in one step
#' hm <- ComplexHeatmap::Heatmap(m,
#'    col=col_hm,
#'    heatmap_legend_param=list(
#'       color_bar="discrete",
#'       border=TRUE,
#'       at=-4:4),
#'    cell_fun=cell_fun_label(m,
#'       col_hm=col_hm))
#' hm
#' 
#' # the cell label function can be created first
#' cell_fun <- cell_fun_label(m,
#'    outline=TRUE,
#'    cex=1.5,
#'    col_hm=col_hm)
#' hm2 <- ComplexHeatmap::Heatmap(m,
#'    col=col_hm,
#'    cell_fun=cell_fun)
#' hm2
#' 
#' @export
cell_fun_label <- function
(m,
   prefix="",
   suffix="",
   cex=1,
   col_hm,
   outline=FALSE,
   abbrev=FALSE,
   ...)
{
   if (!is.list(m)) {
      m <- list(m)
   }
   if (length(prefix) == 0) {
      prefix <- "";
   }
   if (length(suffix) == 0) {
      suffix <- "";
   }
   
   prefix <- rep(prefix, length.out=length(m));
   suffix <- rep(suffix, length.out=length(m));
   cell_fun_return <- function(j, i, x, y, width, height, fill) {
      cell_value <- rmNA(naValue=0, m[[1]][i, j]);
      cell_color <- col_hm(cell_value);
      cell_label <- "";
      for (k in seq_along(m)) {
         mx <- m[[k]]
         cell_value2 <- jamba::rmNA(naValue=0, mx[i, j]);
         if (abbrev && max(cell_value2) >= 1000) {
            cell_label2 <- jamba::asSize(cell_value2,
               kiloSize=1000,
               unitType="",
               sep="");
         } else {
            cell_label2 <- format(cell_value2,
               big.mark=",",
               trim=TRUE)
         }
         cell_label <- paste(cell_label,
            paste0(prefix[k], cell_label2, suffix[k]),
            sep="\n");
      }
      
      if (outline) {
         grid::grid.rect(x=x, y=y,
            width=width,
            height=height,
            gp=grid::gpar(
               col=jamba::setTextContrastColor(cell_color),
               fill=NA))
      }
      if (cex > 0){
         #fontsize <- (10 * cex) - floor(log10(cell_value)) * 1.3;
         fontsize <- (10 * cex);
         grid::grid.text(cell_label, x, y,
            gp=grid::gpar(fontsize=fontsize,
               fontface=1,
               col=jamba::setTextContrastColor(cell_color)
            ));
      }
   }
   cell_fun_return
}

