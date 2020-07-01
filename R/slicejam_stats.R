
## slicejam_stats.R
## - expression matrix operations

#' Normalize SummarizedExperiment data
#' 
#' Normalize SummarizedExperiment data
#' 
#' This function applies one or more data normalization methods
#' to an input `SummarizedExperiment` object. The normalization is
#' applied to one or more matrix data stored in `assays(se)`,
#' each one is run independently.
#' 
#' 
#' @export
se_normalize <- function
(se,
 method=c("quantile", "jammanorm", "limma_batch_adjust"),
 assay_names=NULL,
 genes=NULL,
 samples=NULL,
 params=NULL,
 output_sep="_",
 override=TRUE,
 verbose=FALSE,
 ...)
{
   assay_names <- intersect(assay_names, names(assays(se)));
   if (length(assay_names) == 0) {
      stop("assay_names must be supplied.");
   }
   method <- match.arg(method,
      several.ok=TRUE);
   
   if (length(genes) == 0) {
      genes <- rownames(se);
   } else {
      genes <- intersect(genes, rownames(se));
   }
   if (length(samples) == 0) {
      samples <- colnames(se);
   } else {
      samples <- intersect(samples, colnames(se));
   }
   if (length(genes) == 0 || length(samples) == 0) {
      stop(paste0("Only recognized ",
         length(genes),
         " genes, and ",
         length(samples),
         " samples in the input se."));
   }
   for (assay_name in assay_names) {
      for (imethod in method) {
         output_assay_name <- paste0(
            imethod,
            output_sep,
            assay_name);
         if (output_assay_name %in% names(assays(se)) && !override) {
            if (verbose) {
               jamba::printDebug("se_normalize(): ",
                  sep="",
                  c("Skipped, normalized data already exists for '",
                     output_assay_name,
                     "'"));
            }
            next;
         }
         if (verbose) {
            jamba::printDebug("se_normalize(): ",
               c("Applying method for '",
                  output_assay_name,
                  "'"),
            sep="");
         }
         imatrix <- assays(se[genes, samples])[[assay_name]];
         inorm <- matrix_normalize(imatrix,
            method=imethod,
            params=params);
         assays(se)[[output_assay_name]] <- assays(se)[[assay_name]];
         assays(se)[[output_assay_name]][] <- NA;
         assays(se[genes, samples])[[output_assay_name]] <- inorm[genes, samples];
      }
   }
   return(se);
}

#' Normalize numeric data matrix
#' 
#' Normalize a numeric data matrix
#' 
#' @export
matrix_normalize <- function
(x,
 method=c("quantile", "jammanorm", "limma_batch_adjust"),
 apply_log2=c("ifneeded", "no", "always"),
 floor=0,
 params=list(
   `vsn`=list(lts.quantile=0.5),
   `rsn`=list(excludeFold=2,
      span=0.03),
   cyclicLoess=list(method="default",
      span=0.8),
   `median`=list(Qrange=c(0.5,1.0),
      Irange=NULL,
      HKgenes=NULL,
      useMean=FALSE),
   `jammanorm`=list(controlGenes=NULL,
      minimum_mean=0,
      controlSamples=NULL,
      centerGroups=NULL,
      useMean=TRUE,
      noise_floor=NULL,
      noise_floor_value=NULL),
   `HK`=list(Qrange=NULL,
      Irange=NULL,
      useMean=FALSE,
      HKgenes=c("ACTB", "GAPD", "PPIA", "RPL19")),
   `limma_batch_adjust`=list(
      batch=NULL,
      group=NULL)),
 verbose=FALSE,
 ...)
{
   method <- match.arg(method);
   apply_log2 <- match.arg(apply_log2);

   ## Update the default methodParams without losing the original defaults if not overriden
   params <- update_function_params("matrix_normalize",
      "params",
      params,
      verbose=FALSE);

   if (length(floor) > 0 & any(x < floor)) {
      if (verbose) {
         jamba::printDebug("matrix_normalize(): ",
            c("Applying floor:", floor),
            sep="");
      }
      x[x < floor] <- floor;
   }
   if ("ifneeded" %in% apply_log2) {
      if (any(abs(x) > 40)) {
         x <- jamba::log2signed(x,
            offset=1);
      }
   }
   
   if ("quantile" %in% method) {
      ties <- params$quantile$ties;
      if (length(ties) == 0) {
         ties <- TRUE;
      }
      inorm <- limma::normalizeQuantiles(A=x,
         ties=ties);
   } else if ("jammanorm" %in% method) {
      controlGenes <- params$jammanorm$controlGenes;
      minimum_mean <- params$jammanorm$minimum_mean;
      controlSamples <- params$jammanorm$controlSamples;
      centerGroups <- params$jammanorm$centerGroups;
      useMean <- params$jammanorm$useMean;
      noise_floor <- params$jammanorm$noise_floor;
      noise_floor_value <- params$jammanorm$noise_floor_value;
      if (length(noise_floor) == 0) {
         noise_floor <- -Inf;
      }
      if (length(noise_floor_value) == 0) {
         noise_floor_value <- noise_floor;
      }
      inorm <- jamma::jammanorm(x,
         controlGenes=controlGenes,
         minimum_mean=minimum_mean,
         controlSamples=controlSamples,
         centerGroups=centerGroups,
         useMean=useMean,
         noise_floor=noise_floor,
         noise_floor_value=noise_floor_value,
         verbose=verbose);
   } else if ("limma_batch_adjust" %in% method) {
      batch <- params$limma_batch_adjust$batch;
      group <- params$limma_batch_adjust$group;

      designBatchGroup <- model.matrix(~0+group);
      rownames(designBatchGroup) <- colnames(x);
      
      ## limma::removeBatchEffect()
      inorm <- limma::removeBatchEffect(x,
         batch=batch,
         design=designBatchGroup);
      
   }
   return(inorm);
}

#' Update function default parameters
#' 
#' Update function default parameters
#' 
#' This function is a minor extension to `update_list_elements()`
#' intended to help update function parameters which are defined
#' as a nested list.
#' 
#' @export
update_function_params <- function
(function_name=NULL,
 param_name=NULL,
 new_values=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to facilitate updating default parameters which are present in a list,
   ## but where defaults are defined in a function name.  So if someone wants to change
   ## one of the default values, but keep the rest of the list of defaults, this function
   ## does it.
   ##
   ## The default values are taken from the function formals, using eval(formals(functionName)).
   default_params <- eval(formals(function_name)[[param_name]]);
   if (verbose) {
      jamba::printDebug("update_function_params(): ",
         "str(default_params)");
      print(str(default_params));
   }
   default_params <- update_list_elements(default_params,
      new_values,
      verbose=verbose);
   return(default_params);
}

#' Update a subset of list elements
#' 
#' Update a subset of list elements
#' 
#' This function is intended to help update a nested `source_list`,
#' a subset of whose values should be replaced with entries
#' in `update_list`, leaving any original entries in `source_list`
#' which were not defined in `update_list`.
#' 
#' This function may be useful when manipulating lattice or ggplot2
#' graphical parameters, which are often stored in a nested
#' list structure.
#' 
#' @export
update_list_elements <- function
(source_list,
 update_list,
 list_layer_num=1,
 verbose=TRUE,
 ...)
{
   ## Purpose is to update elements in a list, allowing for multi-layered lists
   ## of lists. In case of a list-of-list, it will call this function again with
   ## each successive layer of listedness.
   ##
   ## Handy for updating lattice graphics settings, which are impossibly nested
   ## tangled ball of textual yarn.  An example:
   ## tp2 <- updateListElements(trellis.par.get(), list(fontsize=list(points=6)));
   ## trellis.par.set(tp2);
   ##
   ## Or in one line:
   ## trellis.par.set(updateListElements(trellis.par.get(), list(fontsize=list(points=6))));
   if (class(update_list) %in% c("list") && class(update_list[[1]]) %in% c("list")) {
      for (update_list_name in names(update_list)) {
         if (update_list_name %in% names(source_list)) {
            ## If the name already exists, we must update items within the list
            source_list[[update_list_name]] <- update_list_elements(
               source_list=source_list[[update_list_name]],
               update_list=update_list[[update_list_name]],
               list_layer_num=list_layer_num+1);
         } else {
            ## If the name does not already exist, we can simply add it.
            source_list[[update_list_name]] <- update_list[[update_list_name]];
         }
      }
   } else {
      if (!is.null(names(update_list))) {
         source_list[names(update_list)] <- update_list;
      } else {
         source_list <- update_list;
      }
   }
   return(source_list);
}

#' Compute contrast statistics on SummarizedExperiment data
#' 
#' Compute contrast statistics on SummarizedExperiment data
#' 
#' This function is essentially a wrapper around statistical methods
#' in the `limma` package, with additional steps to apply statistical
#' thresholds to define "statistical hits" by three main criteria:
#' 
#' * P-value or adjusted P-value
#' * fold change
#' * max group mean
#' 
#' This function is unique in that it applies the statistical methods
#' to one or more "signals" in the input `SummarizedExperiment` assays,
#' specifically intended to compare things like normalization methods.
#' 
#' If multiple statistical thresholds are defined, each one is applied
#' in order, which is specifically designed to compare the effect of
#' applying different statistical thresholds. For example one may want
#' to pre-compute "statistical hits" using adjusted P-value 0.05, and 0.01;
#' or using fold change >= 1.5, or fold change >= 2.0. The underlying
#' statistics are the same, but a column indicating hits is created
#' for each set of thresholds.
#' 
#' Hits are annotated:
#' 
#'  * `-1` for down-regulation
#'  * `0` for un-changed (by the given criteria)
#'  * `1` for up-regulation
#'  
#'  The results are therefore intended to feed directional Venn
#'  diagrams, which display the overlapping gene hits, and whether the
#'  directions are shared or opposed.
#' 
#' This function can optionally apply the limma-voom workflow,
#' which involves calculating matrix weights using `limma::voom()`,
#' then applying those weights during the model fit.
#' 
#' The output is intended to include several convenient formats:
#' 
#' * `statsDFs` - list of `data.frame` stats one per contrast
#' * `statsDF` - one `data.frame` with all stats together
#' * `statsHitsA` - `array` with three dimensions: signal; contrast; threshold
#' whose cells contain hit flags (-1, 0, 1) named by `rownames(se)`.
#' 
#' Design and contrast matrices can be defined using the function
#' `splicejam::groups2contrasts()`. That function assigns each sample
#' to a sample group, then assembles all relevant group contrasts
#' which involve only one-factor contrast at a time. It optionally
#' defines two-factor contrasts (contrast of contrasts) where
#' applicable.
#' 
#' A subset of genes (`rownames(se)`) or samples (`colnames(se)`) can
#' be defined, to restrict calculations to use only the subset data.
#' 
#' @param se `SummarizedExperiment` object
#' @param assay_names `character` vector with one or more assay names
#'    from `names(assays(se))`.
#' 
#' 
#' @export
se_contrast_stats <- function
(se,
   assay_names,
   p_cutoff=0.01,
   adjp_cutoff=NULL,
   fold_cutoff=1.5,
   cutoffAveExpr=NULL,
   mgm_cutoff=NULL,
   p_int_cutoff=cutoffP,
   adjp_int_cutoff=cutoffAdjP,
   confint=FALSE,
   floor_min=NULL,
   floor_value=NULL,
   icontrasts=NULL,
   idesign=NULL,
   igenes=NULL,
   isamples=NULL,
   enforce_design=TRUE,
   use_voom=FALSE,
   weights=NULL,
   robust=FALSE,
   handle_na=c("partial", "full", "none"),
   verbose=TRUE,
   ...)
{
   handle_na <- match.arg(handle_na);
   
   ## Validate input parameters
   if (length(isamples) == 0) {
      isamples <- colnames(se);
   }
   isamples <- intersect(isamples, colnames(se));
   if (length(igenes) == 0) {
      igenes <- rownames(se);
   }
   igenes <- intersect(igenes, rownames(se));

   if (length(idesign) == 0) {
      stop("idesign must be defined.");
   }
   isamples <- intersect(isamples, rownames(idesign));
   idesign <- idesign[match(isamples, rownames(idesign)),,drop=FALSE];
   if (length(rownames(idesign)) == 0) {
      stop("rownames(idesign) must contain values matching colnames(se).");
   }
   icontrasts <- icontrasts[match(rownames(icontrasts), colnames(idesign)),,drop=FALSE];
   if (length(rownames(icontrasts)) == 0) {
      stop("rownames(icontrasts) must match values in colnames(idesign).");
   }
   
   ## Iterate each assay_name
   ## Run statistical tests for gene data
   statsHitsDFs1 <- lapply(jamba::nameVector(assay_names), function(signalSet) {
      retVals <- list();
      imatrix <- SummarizedExperiment::assays(se[igenes, isamples])[[signalSet]];
      if (!"none" %in% handle_na && any(is.na(imatrix))) {
         imatrix <- handle_na_values(imatrix,
            idesign=idesign,
            handle_na=handle_na);
      }
      ## Optionally determine voom weights prior to running limma
      if (use_voom) {
         if (verbose) {
            jamba::printDebug("se_contrast_stats(): ",
               "   Determining Voom weight matrix (within probe reps).");
         }
         imatrix_v <- voom_jam((2^imatrix)-1,
            design=idesign,
            normalize.method="none",
            plot=FALSE,
            verbose=verbose,
            ...);
         weights <- imatrix_v$weights;
         rownames(weights) <- rownames(imatrix);
         colnames(weights) <- colnames(imatrix);
         if (verbose) {
            jamba::printDebug("se_contrast_stats(): ",
               "   Determined Voom weight matrix.");
         }
      }

      #######################################################
      ## Optionally convert zero (or less than zero) to NA
      if (length(floor_min) == 1 && !is.na(floor_min) && any(imatrix <= floor_min)) {
         if (verbose) {
            jamba::printDebug("se_contrast_stats(): ",
               c("Applying floor_min:",
                  floor_min,
                  ", replacing with floor_value:",
                  floor_value),
               sep="");
         }
         imatrix[imatrix <= floor_min] <- floor_value;
      }

      ## Run limma
      rlr_result <- run_limma_replicate(imatrix=imatrix,
         idesign=idesign,
         icontrasts=icontrasts,
         weights=weights,
         robust=robust,
         verbose=verbose,
         ...);
      ## rlr_result is a list
      ## - statsDF
      ## - statsDFs
      ## - repFits=list(subFit1,subFit2,subFit3)
      return(rlr_result);
   });
   
   ## Assemble list of statsDF
   statsDF <- lapply(statsHitsDFs1, function(i){
      i$statsDF;
   });
   ret_list <- list(stats_df=statsDF);
   
   ## Assemble list of statsDFs
   statsHitsDFs <- lapply(statsHitsDFs1, function(i){
      i$statsDFs;
   });
   ret_list$stats_dfs <- statsHitsDFs;
   
   ## list of named lists
   statsHits <- lapply(statsHitsDFs, function(iDFs){
      lapply(iDFs, function(iDF){
         iHitCols <- jamba::nameVector(
            jamba::provigrep("^hit ", colnames(iDF)));
         lapply(iHitCols, function(iHitCol){
            iHitRows <- (!is.na(iDF[,iHitCol]) & iDF[,iHitCol] != 0);
            jamba::nameVector(
               iDF[iHitRows,iHitCol],
               rownames(iDF)[iHitRows]);
         });
      });
   });

   ## array of named lists
   arrayDim <- c(length(statsHits[[1]][[1]]),
      length(statsHits[[1]]),
      length(statsHits));
   arrayDimnames <- list(gsub("[ ]+[^ ]+$", "", names(statsHits[[1]][[1]])), 
      names(statsHits[[1]]), 
      names(statsHits));
   names(arrayDimnames) <- c("HitFilters",
      "Contrasts",
      "Signal");
   jamba::printDebug("arrayDim:", arrayDim);
   jamba::printDebug("arrayDimnames:");print(arrayDimnames);
   statsHitsA <- array(dim=arrayDim,
      data=(unlist(recursive=FALSE,
         unlist(recursive=FALSE,
            statsHits))),
      dimnames=arrayDimnames);
   ret_list$hit_array <- statsHitsA;
   ret_list$hit_list <- statsHits;

   ## Add design and contrast data used
   ret_list$idesign <- idesign;
   ret_list$icontrasts <- icontrasts;

   return(ret_list);   
}

#' Handle NA values in a numeric matrix
#' 
#' Handle NA values in a numeric matrix
#' 
#' @export
handle_na_values <- function
(x,
 idesign,
 handle_na=c("partial", "full", "none"),
 na_value=0,
 na_weight=0,
 verbose=TRUE,
 ...)
{
   handle_na <- match.arg(handle_na);
   xNA <- is.na(x);
   groupL <- multienrichjam::im2list(idesign);
   groupV <- jamba::nameVector(
      rep(names(groupL),
         lengths(groupL)),
      unlist(groupL));
   if ("partial" %in% handle_na) {
      ## We replace NA with zero, except when an entire
      ## group is NA, then we leave it as NA
      # x[xNA] <- 0;
      if (verbose) {
         jamba::printDebug("handle_na_values(): ",
            c("Filling singlet NA with ",
               NAvalue,
               ", leaving full group NA as-is."),
            sep="");
      }
      x <- jamba::rowGroupMeans(x[,names(groupV),drop=FALSE],
         groups=rep(names(groupL), lengths(groupL)),
         rowStatsFunc=function(x,...){
            x1 <- x;
            x1[is.na(x)] <- na_value;
            x1[rowMins(is.na(x)*1) == 1,] <- NA;
            x1;
         });
   } else if ("full" %in% handle_na) {
      ## We replace NA with zero only when an entire group is NA
      ## in order to retain the group in final results.
      ## Otherwise, singlet NA values are kept NA.
      if (verbose) {
         jamba::printDebug("handle_na_values(): ",
            c("Leaving singlet NA as-is, filling full-group NA with ",
               na_value),
            sep="");
      }
      x <- jamba::rowGroupMeans(x[,names(groupV),drop=FALSE],
         groups=rep(names(groupL), lengths(groupL)),
         rowStatsFunc=function(x,...){
            x1 <- x;
            x1[is.na(x)] <- NA;
            x1[rowMins(is.na(x)*1) == 1,] <- na_value;
            x1;
         });
   } else if ("all" %in% handle_na) {
      if (verbose) {
         jamba::printDebug("handle_na_values(): ",
            c("Replacing all NA values with ",
               na_value),
            sep="");
      }
      x[xNA] <- na_value;
   }
   ## define weight matrix
   weights <- jamba::noiseFloor(1-(xNA),
      minimum=na_weight);
   
   return(list(x=x, weights=weights));
}

#' Limma-voom customized for Jam
#' 
#' Limma-voom customized for Jam
#' 
#' This function is based directly upon `limma::voom()` with a
#' small adjustment to handle the presence of `NA` values, which
#' otherwise causes the `stats::lowess()` output to be clearly
#' incorrect. The correction removes `NA` values during this step,
#' producing a result as expected.
#' 
#' @export
voom_jam <- function
(counts,
 design=NULL,
 lib.size=NULL,
 normalize.method="none",
 block=NULL,
 correlation=NULL,
 weights=NULL,
 span=0.5,
 plot=FALSE,
 save.plot=TRUE,
 verbose=FALSE,
 ...)
{
   out <- list()
   
   ## Check counts
   if(is(counts,"DGEList")) {
      out$genes <- counts$genes
      out$targets <- counts$samples
      if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) {
         design <- model.matrix(~group,data=counts$samples)
      }
      if(is.null(lib.size)) {
         lib.size <- with(counts$samples,lib.size*norm.factors);
      }
      counts <- counts$counts;
   } else {
      isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
      if(isExpressionSet) {
         if(length(Biobase::fData(counts))) {
            out$genes <- Biobase::fData(counts)
         }
         if(length(Biobase::pData(counts))) {
            out$targets <- Biobase::pData(counts)
         }
         counts <- Biobase::exprs(counts)
      } else {
         counts <- as.matrix(counts)
      }
   }
   
   n <- nrow(counts);
   if (n < 2L) {
      stop("Need at least two rows to fit a mean-variance trend.");
   }
   
   ## Check design
   if (is.null(design)) {
      design <- matrix(1, ncol(counts), 1);
      rownames(design) <- colnames(counts);
      colnames(design) <- "GrandMean";
   }
   
   ## Check lib.size
   if (is.null(lib.size)) {
      lib.size <- colSums(counts,
         na.rm=TRUE);
   }
   if (verbose) {
      jamba::printDebug("voom_jam(): ",
         "lib.size:", format(digits=2,
            big.mark=",",
            scientific=FALSE,
            lib.size));
      jamba::printDebug("voom_jam(): ",
         "span:", format(digits=2,
            big.mark=",",
            span));
   }
   
   ## Fit linear model to log2-counts-per-million
   y <- t(log2(t(counts + 0.5) / (lib.size + 1) * 1e6));
   y <- limma::normalizeBetweenArrays(y,
      method=normalize.method);
   fit <- limma::lmFit(y,
      design,
      block=block,
      correlation=correlation,
      weights=weights,
      ...);
   if (is.null(fit$Amean)) {
      fit$Amean <- rowMeans(y,
         na.rm=TRUE);
   }
   
   NWithReps <- sum(fit$df.residual > 0L)
   if (NWithReps < 2L) {
      if (NWithReps == 0L) {
         warning("The experimental design has no replication. Setting weights to 1.")
      }
      if (NWithReps == 1L) {
         warning("Only one gene with any replication. Setting weights to 1.")
      }
      out$E <- y;
      out$weights <- y;
      out$weights[] <- 1;
      out$design <- design;
      if (is.null(out$targets)) {
         out$targets <- data.frame(lib.size=lib.size);
      } else {
         out$targets$lib.size <- lib.size;
      }
      return(new("EList", out));
   }
   
   ## Fit lowess trend to sqrt-standard-deviations by log-count-size
   sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e6);
   sy <- sqrt(fit$sigma);
   allzero <- (rowSums(counts, na.rm=TRUE)==0 %in% c(TRUE,NA));
   if (verbose) {
      jamba::printDebug("voom_jam(): ",
         "head(allzero, 10):");
      print(head(allzero, 10));
   }
   if (any(allzero)) {
      sx <- sx[!allzero];
      sy <- sy[!allzero];
   }
   ## 03dec2019: Major change in lowess logic (below)
   ## If there are any NA values in sx or sy, the results
   ## are highly variable, and not correct.
   ## y1 <- rnorm(5000);
   ## y1[sample(1:5000, size=20)] <- NA;
   ## x1 <- seq(from=0.01, to=10, length.out=5000);
   ## plot(x1, y1);
   ## lines(lowess(x1, y1, f=2/3), col="red");
   ## lines(lowess(x1[!is.na(y1)], y1[!is.na(y1)], f=2/3), col="green");
   #l <- stats::lowess(sx,
   #   sy,
   #   f=span);
   s_no_na <- (!is.na(sx) & !is.na(sy));
   l <- stats::lowess(sx[s_no_na],
      sy[s_no_na],
      f=span,
      ...);
   out$sx <- sx;
   out$sy <- sy;
   out$l <- l;
   if (plot) {
      jamba::plotSmoothScatter(x=sx,
         y=sy,
         xlab="log2(count size + 0.5)",
         ylab="sqrt(standard deviation)",
         pch=16,
         cex=0.25);
      title("voom: Mean-variance trend");
      lines(l, col="red");
   }
   
   ## Make interpolating rule
   f <- tryCatch({
      approxfun(l,
         rule=2,
         ties=list("ordered", mean));
   }, error=function(e){
      jamba::printDebug("Error in voom_jam() approxfun():",
         fgText="red");
      #print(e);
      jamba::printDebug("head(fit$Amean, 20):");
      print(head(fit$Amean, 20));
      jamba::printDebug("head(sx, 20):");
      print(head(sx, 20));
      jamba::printDebug("head(fit$sigma, 20):");
      print(head(fit$sigma, 20));
      jamba::printDebug("head(sy, 20):");
      print(head(sy, 20 ));
      jamba::printDebug("str(l):");
      print(str(l));
      stop("voom_jam() failed.");
   });
   
   ## Find individual quarter-root fitted counts
   if (fit$rank < ncol(design)) {
      if (verbose) {
         printDebug("voomJam(): ",
            "Using subset of fit$rank, length(fit$rank):",
            length(fit$rank));
      }
      j <- fit$pivot[seq_len(fit$rank)];
      fitted.values <- fit$coef[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
   } else {
      fitted.values <- fit$coef %*% t(fit$design)
   }
   
   fitted.cpm <- 2^fitted.values;
   fitted.count <- 1e-6 * t(t(fitted.cpm) * (lib.size + 1));
   fitted.logcount <- log2(fitted.count);

   ## Apply trend to individual observations
   w <- 1 / f(fitted.logcount)^4;
   dim(w) <- dim(fitted.logcount);

   ## Output
   rownames(w) <- rownames(y);
   colnames(w) <- colnames(y);
   out$E <- y;
   out$weights <- w;
   out$design <- design;
   out$other <- list(fitted.values=fitted.values);
   if (is.null(out$targets)) {
      out$targets <- data.frame(lib.size=lib.size);
   } else {
      out$targets$lib.size <- lib.size;
   }
   if (save.plot) {
      out$voom.xy <- list(x=sx,
         y=sy,
         xlab="log2( count size + 0.5 )",
         ylab="Sqrt( standard deviation )");
      out$voom.line <- l;
   }
   new("EList", out);
}

#' Run limma contrasts with optional probe replicates
#' 
#' @export
run_limma_replicate <- function
(imatrix,
 idesign,
 icontrasts,
 weights=NULL,
 robust=FALSE,
 adjust.method="BH",
 confint=FALSE,
 verbose=FALSE,
 ...)
{
   imatrixES <- Biobase::ExpressionSet(assayData=imatrix,
      featureData=new("AnnotatedDataFrame",
         data=data.frame(
            probes=rownames(imatrix),
            row.names=rownames(imatrix)
         )
      )
   );

   ## lmFit
   subFit1 <- limma::lmFit(imatrixES,
      design=idesign,
      weights=weights);
   if (length(rownames(subFit1$coefficients)) == 0) {
      rownames(subFit1$coefficients) <- rownames(subFit1$genes);
   }

   ## Add back in the Amean values (probe means)
   if (length(subFit1$Amean) == 0) {
      subFit1$Amean <- jamba::nameVector(
         rowMeans(imatrix, na.rm=TRUE),
         rownames(imatrix));
      if (verbose) {
         jamba::printDebug("run_limma_replicate(): ",
            "Added Amean values since lmFit did not.");
      }
   }

   ## run contrasts on the limma model
   if (verbose) {
      jamba::printDebug("run_limma_replicate(): ",
         "Running contrasts.fit on subFit1");
   }
   subFit2 <- limma::contrasts.fit(subFit1,
      icontrasts);

   ## eBayes adjustment for signal-based noise
   if (verbose) {
      jamba::printDebug("run_limma_replicate(): ",
         "Running eBayes on subFit2");
   }
   subFit3 <- limma::eBayes(subFit2,
      robust=robust);
   #if (nrow(subMatrix) == ndups) {
   #   subFit3 <- subFit3[1,];
   #}
   
   ## define colnames to rename to include the contrast name
   renameCols <- c("logFC",
      "P.Value",
      "adj.P.Val",
      "t",
      "F",
      "B",
      "CI.L",
      "CI.R",
      "AveExpr");
   
   ## Get summary table for each contrast
   contrastNames <- colnames(subFit3$contrast);
   
   dimNum <- nrow(subFit3$coefficients);
   ## top table for each contrast
   tts <- lapply(jamba::nameVector(contrastNames), function(contrastName) {
      tt <- limma::topTable(subFit3,
         coef=contrastName,
         number=Inf,
         adjust.method=adjust.method,
         confint=confint);
      
      ## Rename columns to include the contrast
      ## Note: this step could rename "probes" or "genes"
      tt <- jamba::renameColumn(tt,
         from=c(renameCols),
         to=c(paste(renameCols, contrastName)));
      
      if (verbose) {
         jamba::printDebug("run_limma_replicate(): ",
            "contrastName:",
            contrastName,
            " head(topTable):");
         print(head(tt, 5));
      }
      tt;
   } );
   
   # Create stats table overall
   ttAll1 <- topTable(subFit3,
      number=Inf,
      adjust.method=adjust.method,
      confint=confint);
   ttAll <- ttAll1[match(rownames(imatrix), ttAll1[,1]),,drop=FALSE];
   ## 30jun2020: changed to keep empty rows
   ttAll[,1] <- rownames(imatrix);
   #ttAll <- ttAll[!is.na(ttAll[,1]),,drop=FALSE];

   ## Rename columns to include the contrast
   ttAll <- jamba::renameColumn(ttAll,
      from=c(renameCols),
      to=c(paste(renameCols)));

   #tts <- list(overall=ttAll,
   #   list=tts);
   #all_fits <- list(tts=tts,
   #   subFit3=subFit3,
   #   subFit1=subFit1);

   trim_colnames <- c("t", "B", "F");
   trim_grep <- paste0("^", trim_colnames, "( |[(]|$)");
   trim_grep;
   statsDFs <- lapply(contrastNames, function(contrastName){
      statsDF <- tts[[contrastName]];
      trimCols1 <- jamba::provigrep(trim_grep, colnames(statsDF));
      if (length(trimCols1) > 0) {
         statsDF <- statsDF[,setdiff(colnames(statsDF), trimCols1),drop=FALSE];
      }
      if (length(jamba::tcount(statsDF[,1], minCount=2)) == 0) {
         rownames(statsDF) <- statsDF[,1];
      }
      statsDF;
   });
   
   all_statsDF_colnames <- unique(c(colnames(ttAll),
      unlist(lapply(statsDFs, colnames))));
   all_statsDF <- jamba::mergeAllXY(c(
      list(overall=ttAll),
      statsDFs));
   all_statsDF <- all_statsDF[, all_statsDF_colnames, drop=FALSE];
   if (length(jamba::tcount(all_statsDF[,1], minCount=2)) == 0) {
      rownames(all_statsDF) <- all_statsDF[,1];
   }
   return(
      list(statsDF=all_statsDF,
         statsDFs=statsDFs,
         repFits=list(subFit1=subFit1,
            subFit2=subFit2,
            subFit3=subFit3)
      )
   )
}
