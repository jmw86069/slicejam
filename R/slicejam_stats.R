
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
#' Note that supplying `genes` and `samples` will apply normalization
#' to only those `genes` and `samples`, and this data will be
#' stored in the full `SummarizedExperiment` object `se` with
#' `NA` values used to fill any values not present in `genes`
#' or `samples`.
#' 
#' For example if `assay_names` contains two assay names,
#' and `method` contains two methods, the output will include
#' four normalizations, where each assay name is normalized two ways.
#' The output assay names will be something like `"assay1_method1"`,
#' `"assay1_method2"`, `"assay2_method1"`, `"assay2_method2"`.
#' It is not always necessary to normalize data by multiple different
#' methods, however when two methods are similar and need to be
#' compared, the `SummarizedExperiment` object is a convenient
#' place to store different normalization results for downstream
#' comparison. Further, the method `se_contrast_stats()` is able
#' to apply equivalent statistical contrasts to each normalization,
#' and returns an array of statistical hits which is convenient
#' for direct comparison of results.
#' 
#' This method calls `matrix_normalize()` to perform each normalization
#' step, see that function description for details on each method.
#' 
#' @family slicejam stats
#' 
#' @return `SummarizedExperiment` object where the normalized output
#'    is added to `assays(se)` using the naming format `method_assayname`.
#' 
#' @param se `SummarizedExperiment` object
#' @param method `character` vector indicating which normalization method(s)
#'    to apply.
#' @param assay_names `character` vector or one or more `names(assays(se))`
#'    that indicates which numeric matrix to use during normalization. When
#'    multiple values are provided, each matrix is normalized independently
#'    by each `method`.
#' @param genes `character` vector (optional) used to define a subset of
#'    gene rows in `se` to use for normalization.
#'    Values must match `rownames(se)`.
#' @param samples `character vector (optional) used to define a subset of
#'    sample columns in `se` to use for normalization.
#'    Values must match `colnames(se)`.
#' @param params `list` (optional) parameters specific to each
#'    normalization method, passed to `matrix_normalize()`. Any 
#'    value which is not defined in the `params` provided will use
#'    the default value in `matrix_normalize()`, for example
#'    `params=list(jammanorm=list(minimum_mean=2))` will use
#'    `minimum_mean=2` then use other default values relevant
#'    to the `jammanorm` normalization method.
#' @param output_sep `character` string used as a delimited between the
#'    `method` and the `assay_names` to define the output assay name,
#'    for example when `assay_name="counts"`, `method="quantile"`,
#'    and `output_sep="_"` the new assay name will be `"quantile_counts"`.
#' @param override `logical` indicating whether to override any pre-existing
#'    matrix values with the same output assay name. When `override=FALSE`
#'    and the output assay name already exists, the normalization will
#'    not be performed.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `matrix_normalize()`.
#' 
#' @examples
#' if (jamba::check_pkg_installed("farrisdata")) {
#' 
#'    # se_normalize
#'    suppressPackageStartupMessages(library(SummarizedExperiment))
#'    GeneSE <- farrisdata::farrisGeneSE;
#'    samples <- colnames(GeneSE);
#'    genes <- rownames(GeneSE);
#'    
#'    GeneSE <- se_normalize(GeneSE,
#'       genes=genes,
#'       samples=samples,
#'       assay_names=c("raw_counts", "counts"),
#'       method="jammanorm",
#'       params=list(jammanorm=list(minimum_mean=5)))
#'    names(assays(GeneSE))
#'    
#'    # review normalization factor values
#'    round(digits=3, attr(assays(GeneSE)$jammanorm_raw_counts, "nf"))
#'    
#'    # the data in "counts" was already normalized
#'    # so the normalization factors are very near 0 as expected
#'    round(digits=3, attr(assays(GeneSE)$jammanorm_counts, "nf"))
#'    
#'    
#'    # note that housekeeper genes are supplied in params
#'    set.seed(123);
#'    hkgenes <- sample(rownames(GeneSE), 1000)
#'    GeneSE <- se_normalize(GeneSE,
#'       genes=genes,
#'       samples=samples,
#'       assay_names=c("raw_counts"),
#'       method="jammanorm",
#'       params=list(jammanorm=list(minimum_mean=5,
#'          controlGenes=hkgenes)))
#'    round(digits=3, attr(assays(GeneSE)$jammanorm_raw_counts, "nf"))
#'    
#'    # example showing quantile normalization
#'    GeneSE <- se_normalize(GeneSE,
#'       assay_names=c("raw_counts"),
#'       method="quantile",
#'       params=list(jammanorm=list(min_mean=5)))
#' }
#'    
#' 
#' @export
se_normalize <- function
(se,
 method=c("quantile",
    "jammanorm",
    "limma_batch_adjust"),
 assay_names=NULL,
 genes=NULL,
 samples=NULL,
 params=list(
 `quantile`=list(
     ties=TRUE),
 `jammanorm`=list(controlGenes=NULL,
     minimum_mean=0,
     controlSamples=NULL,
     centerGroups=NULL,
     useMedian=FALSE,
     noise_floor=NULL,
     noise_floor_value=NULL),
 `limma_batch_adjust`=list(
     batch=NULL,
     group=NULL)),
 output_sep="_",
 override=TRUE,
 verbose=FALSE,
 ...)
{
   assay_names <- intersect(assay_names, names(SummarizedExperiment::assays(se)));
   if (length(assay_names) == 0) {
      stop("assay_names must be supplied.");
   }
   method <- match.arg(method,
      several.ok=TRUE);
   
   allgenes <- rownames(se);
   allsamples <- colnames(se);
   if (length(genes) == 0) {
      genes <- allgenes;
   } else {
      genes <- intersect(genes, allgenes);
   }
   if (length(samples) == 0) {
      samples <- allsamples;
   } else {
      samples <- intersect(samples, allsamples);
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
         if (output_assay_name %in% names(SummarizedExperiment::assays(se)) && !override) {
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
         imatrix <- SummarizedExperiment::assays(se[genes, samples])[[assay_name]];
         inorm <- matrix_normalize(imatrix,
            method=imethod,
            params=params,
            verbose=verbose - 1,
            ...);
         
         # generate matrix of NA values to fill for normalized genes, samples
         # this way we can add attributes to the full matrix
         # when normalizing a smaller matrix
         namatrix <- matrix(data=NA,
            ncol=length(allsamples),
            nrow=length(allgenes),
            dimnames=list(allgenes,
               allsamples));
         # assign normalized data for genes, samples
         namatrix[genes, samples] <- inorm;
         
         # re-assign any missing attributes
         # potentially useful information from the normalization method
         na_attrnames <- names(attributes(namatrix));
         new_attrs <- setdiff(names(attributes(inorm)), na_attrnames);
         if (length(new_attrs) > 0) {
            if (verbose) {
               jamba::printDebug("   re-assign new_attrs: ",
                  sep=", ",
                  new_attrs);
            }
            for (new_attr in new_attrs) {
               attr(namatrix, new_attr) <- attr(inorm, new_attr);
            }
         }
         
         # assign namatrix to se
         SummarizedExperiment::assays(se)[[output_assay_name]] <- namatrix;
      }
   }
   return(se);
}

#' Normalize a numeric data matrix
#' 
#' Normalize a numeric data matrix
#' 
#' This function is a wrapper for several relevant normalization
#' methods that operate on a numeric matrix.
#' 
#' # Normalization Methods Implemented:
#' 
#' ## method='quantile'
#' 
#' Quantile-normalization performed by
#' `limma::normalizeQuantiles()`. This method has one
#' parameter `"ties"` passed to `limma::normalizeQuantiles()`,
#' the default here `ties=TRUE` which handles tied numeric
#' expression values in a robust way to avoid unpredictability
#' otherwise. This option is especially relevant with expression
#' count data, where integer counts cause a large number
#' of values to be represented multiple times.
#' 
#' ## method='jammanorm'
#' 
#' Median-normalization performed by
#' `jamma::jammanorm()`. This method shifts expression
#' data as shown on MA-plots, so the median expression
#' is zero across all samples, using only the rows that
#' meet the relevant criteria.
#' 
#' Some relevant criteria to
#' define rows used for normalization:
#' 
#' * `controlGenes` defines specific genes to use for
#' normalization, such as housekeeper genes. It may also
#' be useful to use detected genes here, so the normalization
#' only considers those genes defined as detected by
#' the protocol.
#' * `minimum_mean` sets a numeric threshold and requires
#' the mean expression (shown on the x-axis of the MA-plot)
#' to be at least this value.
#' 
#' Note that when both `controlGenes` and `minimum_mean`
#' are defined, both criteria are enforced. So the `controlGenes`
#' are also required to have expression of at least `minimum_mean`.
#' 
#' Also note that all rows of data are normalized by this method,
#' only the subset of rows defined by `controlGenes` and `minimum_mean`
#' are used to compute the normalization factor.
#' 
#' 
#' ## method='limma_batch_adjust'
#' 
#' Batch adjustment performed by
#' `limma::removeBatchEffect()` which is intended to apply
#' batch-adjustment as a form of normalization, but which
#' does not represent full normalization itself. There are
#' two relevant parameters: `"batch"` which is a vector of
#' batch values in order of `colnames(x)`, and `"group"`
#' which is a vector of sample groups in order of `colnames(x)`.
#' 
#' # Other useful parameters
#' 
#' Note the `floor` and `enforce_norm_floor` have recommended
#' default values `floor=0` and `enforce_norm_floor=TRUE`. These
#' defaults will set any assay value at or below `0` to `0`,
#' and after normalization any values whose input values were
#' at or below `0` will also be set to `0` to prevent normalizing
#' a value of `0` to non-zero. Any normalized value at or
#' below `0` will also be set to `0` to prevent results from
#' containing negative normalized values.
#' 
#' The assumption for this default is that a value of zero
#' is not a measurement but represents the lack of a measurement.
#' Similarly, the intent of `floor` is a numeric threshold at or
#' below there is no confidence in the reported measurement, therefore
#' values at or below this threshold are treated as equivalent
#' to the threshold for the purpose of downstream analyses.
#' 
#' Some platforms like QPCR for example, have substantially lower
#' confidence at high CT values, where expression values
#' using the equation `2^(40-CT)` might impose a noise threshold
#' at expression 32 or lower. This noise threshold for QPCR
#' means any expression measurement of 32 or lower is as likely
#' to be `32` as it is to be `2`, and therefore any differences
#' between reported expression of `32` and `2` should not be
#' considered relevant. Applying `floor=32` in this case
#' accomplishes this goal by setting all values at or below
#' `32` to `32`. Of course when using this method `matrix_normalize()`
#' the data should be log2 transformed, which means the `floor`
#' should also be log2 transformed, e.g. `floor=log2(32)`
#' which is `floor=5`.
#' 
#' One alternative might be to set values at or below zero to `NA`
#' prior to normalization, and before calling `matrix_normalize()`.
#' In this case, only non-NA values will be used during
#' normalization according to the `method` being used.
#' 
#' @return `numeric` matrix with the same dimensions as the
#'    input matrix `x`. Some normalization methods return
#'    additional information in `attributes(x)`, for example
#'    `method="jammanorm"` will return the vector of housekeeper
#'    genes used in `attr(x, "hk")` for normalization of each sample
#'    when supplied with `controlGenes` values.
#' 
#' @family slicejam stats
#' 
#' @param x `numeric` matrix with sample columns, and typically
#'    gene rows, but any measured assay row will meet the assumptions
#'    of the method.
#' @param method `character` string indicating which normalization
#'    method to apply.
#' @param apply_log2 `character` string indicating whether to apply
#'    log2 transformation: `"ifneeded"` will apply log2 transform
#'    when any absolute value is greater than 40; `"no"` will not
#'    apply log2 transformation; `"always"` will apply log2 transform.
#'    Note the log2 transform is applied with `jamba::log2signed(x, offset=1)`
#'    which is equivalent to `log(1 + x)` except that negative values
#'    are also transformed using the absolute value, then multiplied
#'    by their original sign.
#' @param floor `numeric` value indicating the lowest accepted numeric
#'    value, below which values are assigned to this floor. The default
#'    `floor=0` requires all values are `0`, and any values below `0` are
#'    assigned `0`. Note that the `floor` is applied after log2 transform,
#'    when the log2 transform is performed.
#' @param enforce_norm_floor `logical` indicating whether to enforce the
#'    `floor` for the normalized results, default is `TRUE`. For example,
#'    when `floor=0` any values at or below `0` are set to `0` before
#'    normalization. After normalization some of these values will be
#'    above or below `0`. When `enforce_norm_floor=TRUE` these values
#'    will again be set to `0` because they are considered to be
#'    below the noise threshold of the protocol, and adjustments
#'    are not relevant; also any normalized values below the `floor`
#'    will also be set to `floor`.
#' @param params `list` of parameters relevant to the `method` of
#'    normalization. The `params` should be a `list` named by the `method`,
#'    whose values are a list named by the relevant method parameter.
#'    See examples.
#' @param verbose `logical` indicating whether to print verbose output.
#' 
#' @examples
#' # use farrisdata real world data if available
#' if (jamba::check_pkg_installed("farrisdata")) {
#' 
#'    suppressPackageStartupMessages(library(SummarizedExperiment))
#'    
#'    # test matrix_normalize()
#'    GeneSE <- farrisdata::farrisGeneSE;
#'    imatrix <- assays(GeneSE)$raw_counts;
#'    genes <- rownames(imatrix);
#'    samples <- colnames(imatrix);
#'    head(imatrix);
#'    
#'    # matrix_normalize()
#'    # normalize the numeric matrix directly
#'    imatrix_norm <- matrix_normalize(imatrix,
#'       genes=genes,
#'       samples=samples,
#'       method="jammanorm",
#'       params=list(minimum_mean=5))
#'    names(attributes(imatrix_norm))
#'    
#'    # review normalization factors
#'    round(digits=3, attr(imatrix_norm, "nf"));
#'    
#'    # example for quantile normalization
#'    imatrix_quant <- matrix_normalize(imatrix,
#'       genes=genes,
#'       samples=samples,
#'       method="quantile")
#'    names(attributes(imatrix_quant))
#' }
#' 
#' 
#' # simulate reasonably common expression matrix
#' set.seed(123);
#' x <- matrix(rnorm(9000)/4, ncol=9);
#' colnames(x) <- paste0("sample", LETTERS[1:9]);
#' rownames(x) <- paste0("gene", jamba::padInteger(seq_len(nrow(x))))
#' rowmeans <- rbeta(nrow(x), shape1=2, shape2=5)*14+2;
#' x <- x + rowmeans;
#' for (i in 1:9) {
#'    x[,i] <- x[,i] + rnorm(1);
#' }
#' 
#' # display MA-plot with jamma::jammaplot()
#' jamma::jammaplot(x)
#' 
#' # normalize by jammanorm
#' xnorm <- matrix_normalize(x, method="jammanorm")
#' jamma::jammaplot(xnorm, maintitle="method='jammanorm'")
#' 
#' # normalize by jammanorm with housekeeper genes
#' hk_genes <- sample(rownames(x), 10);
#' xnormhk <- matrix_normalize(x,
#'    method="jammanorm",
#'    params=list(jammanorm=list(controlGenes=hk_genes)))
#'    
#' jamma::jammaplot(xnormhk,
#'    maintitle="method='jammanorm' with housekeeper genes",
#'    highlightPoints=list(housekeepers=hk_genes),
#'    highlightColor="green");
#' 
#' xnormhk6 <- matrix_normalize(x,
#'    method="jammanorm",
#'    params=list(jammanorm=list(
#'       controlGenes=hk_genes,
#'       minimum_mean=6)))
#' hk_used <- attr(xnormhk6, "hk")[[1]];
#' jamma::jammaplot(xnormhk6,
#'    maintitle="method='jammanorm' with housekeeper genes, minimum_mean=6",
#'    highlightPoints=list(housekeepers=hk_genes,
#'       hk_used=hk_used),
#'    highlightColor=c("red", "green"));
#' 
#' # normalize by quantile
#' xquant <- matrix_normalize(x, method="quantile")
#' jamma::jammaplot(xquant,
#'    maintitle="method='quantile'")
#' 
#' # simulate higher noise for lower signal
#' rownoise <- rnorm(prod(dim(x))) * (3 / ((rowmeans*1.5) - 1.5));
#' xnoise <- x;
#' xnoise <- xnoise + rownoise;
#' jamma::jammaplot(xnoise,
#'    maintitle="simulated higher noise at lower signal");
#' 
#' # simulate non-linearity across signal
#' # sin(seq(from=pi*4/10, to=pi*7/10, length.out=100))-0.8
#' rowadjust <- (sin(pi * jamba::normScale(rowmeans, from=3.5/10, to=5.5/10)) -0.9) * 20;
#' xwarp <- xnoise;
#' xwarp[,3] <- xnoise[,3] + rowadjust;
#' jamma::jammaplot(xwarp,
#'    maintitle="signal-dependent noise and non-linear effects");
#' 
#' # quantile-normalization is indicated for this scenario
#' xwarpnorm <- matrix_normalize(xwarp,
#'    method="quantile");
#' jp <- jamma::jammaplot(xwarpnorm,
#'    maintitle="quantile-normalized: signal-dependent noise and non-linear effects");
#' 
#' @export
matrix_normalize <- function
(x,
 method=c("quantile", "jammanorm", "limma_batch_adjust"),
 apply_log2=c("ifneeded", "no", "always"),
 floor=0,
 enforce_norm_floor=TRUE,
 params=list(
   #`vsn`=list(lts.quantile=0.5),
   #`rsn`=list(excludeFold=2,
   #   span=0.03),
   #cyclicLoess=list(method="default",
   #   span=0.8),
   `quantile`=list(
      ties=TRUE),
   `jammanorm`=list(controlGenes=NULL,
      minimum_mean=0,
      controlSamples=NULL,
      centerGroups=NULL,
      useMedian=FALSE,
      noise_floor=NULL,
      noise_floor_value=NULL),
   `limma_batch_adjust`=list(
      batch=NULL,
      group=NULL)),
 verbose=TRUE,
 ...)
{
   method <- match.arg(method);
   apply_log2 <- match.arg(apply_log2);

   ## Update the default methodParams without losing the original defaults if not overriden
   params <- update_function_params("matrix_normalize",
      "params",
      params,
      verbose=FALSE);

   # apply log2 transform if needed
   if ("ifneeded" %in% apply_log2) {
      if (any(abs(x) > 40)) {
         x <- jamba::log2signed(x,
            offset=1);
         if (verbose) {
            jamba::printDebug("matrix_normalize(): ",
               c("Applied ",
                  "jamba::log2signed(x, offset=1)",
                  " because ",
                  "any(abs(x) > 40)"),
               sep="");
         }
      }
   } else if ("always" %in% apply_log2) {
      x <- jamba::log2signed(x,
         offset=1);
      if (verbose) {
         jamba::printDebug("matrix_normalize(): ",
            c("Applied ", "jamba::log2signed(x, offset=1)"),
            sep="");
      }
   }
   
   # apply optional floor
   if (length(floor) > 0 & any(x <= floor)) {
      x_floored <- (x <= floor);
      x[x_floored] <- floor;
      if (verbose) {
         jamba::printDebug("matrix_normalize(): ",
            c("Applied floor:", floor),
            sep="");
      }
   } else {
      x_floored <- FALSE
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
      useMedian <- params$jammanorm$useMedian;
      noise_floor <- params$jammanorm$noise_floor;
      noise_floor_value <- params$jammanorm$noise_floor_value;
      if (length(noise_floor) == 0) {
         noise_floor <- -Inf;
      }
      if (length(noise_floor_value) == 0) {
         noise_floor_value <- noise_floor;
      }
      if (verbose) {
         jamba::printDebug("matrix_normalize(): ",
            "Calling ", "jammanorm():",
               c("\n      minimum_mean:", minimum_mean,
               "\n      useMedian:", useMedian,
               if (length(noise_floor) > 0 && noise_floor > -Inf) {
                  c("\n      noise_floor:", noise_floor)},
               if ((length(jamba::rmNA(noise_floor_value)) > 0 && jamba::rmNA(noise_floor_value) > -Inf) ||
                     any(is.na(noise_floor_value))) {
                  c("\n      noise_floor_value:", noise_floor_value)}),
               sep="");
      }
      inorm <- jamma::jammanorm(x,
         controlGenes=controlGenes,
         minimum_mean=minimum_mean,
         controlSamples=controlSamples,
         centerGroups=centerGroups,
         useMedian=useMedian,
         useMean=NULL,
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
   
   # enforce the floor for output matrix
   if (enforce_norm_floor && any(x_floored)) {
      inorm[x_floored | inorm <= floor] <- floor;
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
#' @family slicejam utilities
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
#' @family slicejam utilities
#' 
#' @export
update_list_elements <- function
(source_list,
 update_list,
 list_layer_num=1,
 verbose=FALSE,
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
   if (length(update_list) == 0) {
      return(source_list);
   }
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
#' @family slicejam stats
#' 
#' @export
se_contrast_stats <- function
(se,
 assay_names,
 adjp_cutoff=0.05,
 p_cutoff=NULL,
 fold_cutoff=1.5,
 int_adjp_cutoff=adjp_cutoff,
 int_p_cutoff=p_cutoff,
 int_fold_cutoff=fold_cutoff,
 mgm_cutoff=NULL,
 ave_cutoff=NULL,
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
 collapse_by_gene=FALSE,
 verbose=FALSE,
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
   stats_hits_dfs1 <- lapply(jamba::nameVector(assay_names), function(signalSet) {
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
         adjp_cutoff=adjp_cutoff,
         p_cutoff=p_cutoff,
         fold_cutoff=fold_cutoff,
         mgm_cutoff=mgm_cutoff,
         int_adjp_cutoff=int_adjp_cutoff,
         int_p_cutoff=int_p_cutoff,
         int_fold_cutoff=int_fold_cutoff,
         ave_cutoff=ave_cutoff,
         collapse_by_gene=collapse_by_gene,
         ...);
      ## rlr_result is a list
      ## - statsDF
      ## - statsDFs
      ## - repFits=list(subFit1,subFit2,subFit3)
      return(rlr_result);
   });
   
   ## Assemble list of statsDF
   stats_df <- lapply(stats_hits_dfs1, function(i){
      i$stats_df;
   });
   ret_list <- list(stats_df=stats_df);
   
   ## Assemble list of statsDFs
   stats_dfs <- lapply(stats_hits_dfs1, function(i){
      i$stats_dfs;
   });
   ret_list$stats_dfs <- stats_dfs;
   
   ## list of named lists
   stats_hits <- lapply(stats_dfs, function(iDFs){
      lapply(iDFs, function(iDF){
         iHitCols <- jamba::nameVector(
            jamba::provigrep("^hit", colnames(iDF)));
         if (verbose) {
            jamba::printDebug("se_contrast_stats(): ",
               "iHitCols:", iHitCols);
         }
         lapply(iHitCols, function(iHitCol){
            iHitRows <- (!is.na(iDF[,iHitCol]) & iDF[,iHitCol] != 0);
            jamba::nameVector(
               iDF[iHitRows,iHitCol],
               rownames(iDF)[iHitRows]);
         });
      });
   });
   if (verbose) {
      jamba::printDebug("se_contrast_stats(): ",
         "ssdim(stats_hits[[1]]):");
      print(jamba::ssdim(stats_hits[[1]]));
   }

   ## array of named lists
   arrayDim <- c(length(stats_hits[[1]][[1]]),
      length(stats_hits[[1]]),
      length(stats_hits));
   arrayDimnames <- list(gsub("[ ]+[^ ]+$", "", names(stats_hits[[1]][[1]])), 
      names(stats_hits[[1]]), 
      names(stats_hits));
   names(arrayDimnames) <- c("Cutoffs",
      "Contrasts",
      "Signal");
   if (verbose) {
      jamba::printDebug("se_contrast_stats(): ",
         "arrayDim:",
         arrayDim);
      jamba::printDebug("se_contrast_stats(): ",
         "arrayDimnames:");
      print(arrayDimnames);
   }
   hit_array <- array(dim=arrayDim,
      data=(unlist(recursive=FALSE,
         unlist(recursive=FALSE,
            stats_hits))),
      dimnames=arrayDimnames);
   ret_list$hit_array <- hit_array;
   ret_list$hit_list <- stats_hits;

   ## Add design and contrast data used
   ret_list$idesign <- idesign;
   ret_list$icontrasts <- icontrasts;

   return(ret_list);   
}

#' Handle NA values in a numeric matrix
#' 
#' Handle NA values in a numeric matrix
#' 
#' @family slicejam stats
#' 
#' @export
handle_na_values <- function
(x,
 idesign,
 handle_na=c("partial", "full", "none"),
 na_value=0,
 na_weight=0,
 verbose=FALSE,
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
#' @family slicejam stats
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
         jamba::printDebug("voom_jam(): ",
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
#' Run limma contrasts with optional probe replicates
#' 
#' This function is called by `se_contrast_stats()` to perform
#' the comparisons defined as contrasts. The `se_contrast_stats()`
#' function operates on a `SummarizedExperiment` object,
#' this function operates on the `numeric` `matrix` values
#' directly.
#' 
#' This function also calls `ebayes2dfs()` which extracts
#' each contrast result as a `data.frame`, whose column names
#' are modified to include the contrast names.
#' 
#' This function optionally (not yet ported from previous
#' implementation) detects replicate probes, and performs
#' the internal correlation calculations recommended by
#' `limma user guide` for replicate probes. In that case,
#' it detects each level of probe replication so that
#' each can be properly calculated. For example, Agilent
#' human 4x44 arrays often contain a large number of probes
#' with 8 replicates; a subset of probes with 4 replicates;
#' then the remaining probes (the majority overall) have
#' only one replicate each. In that case, this function
#' splits data into 8-replicate, 4-replicate, and 1-replicate
#' subsets, calculates correlations on the 8-replicate and
#' 4-replicate subsets separately, then runs limma calculations
#' on the three subsets independently, then merges the results
#' into one large table. The end result is that the
#' final table contains one row per unique probe after
#' adjusting for probe replication properly in each scenario.
#' As the Agilent microarray layout is markedly less widely
#' used that in past, the priority to port this methodology
#' is quite low.
#' 
#' 
#' @family slicejam stats
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
 trim_colnames=c("t", "B", "F"),
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
   if (1 == 1) {
      stats_dfs <- ebayes2dfs(lmFit3=subFit3,
         lmFit1=subFit1,
         define_hits=TRUE,
         trim_colnames=trim_colnames,
         verbose=verbose,
         ...);
   } else {
      stats_dfs <- lapply(jamba::nameVector(contrastNames), function(contrastName) {
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
   }
   
   stats_df <- jamba::mergeAllXY(stats_dfs);
   stats_df_colnames <- unique(
      unlist(lapply(stats_dfs, colnames)));
   stats_df <- jamba::mergeAllXY(stats_dfs);
   stats_df <- stats_df[, stats_df_colnames, drop=FALSE];
   
   if (length(jamba::tcount(stats_df[,1], minCount=2)) == 0) {
      rownames(stats_df) <- stats_df[,1];
   }
   return(
      list(stats_df=stats_df,
         stats_dfs=stats_dfs,
         rep_fits=list(lmFit1=subFit1,
            lmFit2=subFit2,
            lmFit3=subFit3)
      )
   )
}


#' Convert limma eBayes fit to data.frame with annotated hits
#' 
#' Convert limma eBayes fit to data.frame with annotated hits
#' 
#' This function is called by `run_limma_replicate()` as
#' an extension to `limma::topTable()`, that differs in that
#' it is performed for each contrast in the input `lmFit3` object.
#' 
#' By default the columns include the contrast, so that each `data.frame`
#' is self-described.
#' 
#' When `define_hits=TRUE`, then statistical thresholds are applied
#' to define a set of statistical hits. The thresholds available include:
#' 
#' 1. `adjp_cutoff` - applied to `"adj.P.Val"` for adjusted P-value.
#' 2. `p_cutoff` - applied to `"P.Value"` for raw, unadjusted P-value.
#' 3. `fold_cutoff` - normal space fold change, applied to `"logFC"`
#'    by using `log2(fold_cutoff)`.
#' 4. `mgm_cutoff` - max group mean, applied to the highest group mean
#'    value involved in each specific contrast.
#' 5. `ave_cutoff` - applied to `"AveExpr"` which represents the mean
#'    value across all sample groups.
#' 
#' Note that `mgm_cutoff` requires input `lmFit1` which stores the
#' group mean values used in the limma workflow.
#' 
#' Note also there are optional arguments specific to interaction
#' contrasts, which in this context is assumed to be a
#' "fold change of fold changes" style of contrast, for example:
#' `(groupA-groupB)-(groupC-groupD)`. The purpose is distinct interaction
#' thresholds is to enable reasonable data mining, sometimes with
#' somewhat more lenient thresholds for interaction contrasts.
#' For example, one may use `adjp_cutoff=0.01` and `int_adjp_cutoff=0.05`,
#' or `fold_cutoff=2` and `int_fold_cutoff=1.5`.
#' 
#' By default, `rename_headers=TRUE` causes colnames to include the
#' contrast, for example renaming colname `"logFC"` to `"logFC contrastA"`.
#' This change helps reinforce the source of the statistical results,
#' and allows the `data.frame` results to be merged together using
#' `base::merge()`.
#' 
#' Indeed, `merge_df=TRUE` will cause all `data.frame` results to be
#' merged into one large `data.frame`, using `jamba::mergeAllXY()`.
#' 
#' @return `list` with one `data.frame` per contrast defined in
#'    the input `lmFit3` object. When `define_hits=TRUE` there
#'    will be one column per statistical threshold, named `"hit"`
#'    followed by an abbreviation of the statistical thresholds
#'    which were applied.
#'    When `merge_df=TRUE` the returned data will be one
#'    `data.frame` object.
#' 
#' @family slicejam stats
#' 
#' @param lmFit3 object returned by `limma::eBayes()`.
#' @param lmFit1 object returned by `limma::lmFit()`, optional.
#' @param define_hits `logical` indicating whether to define hits
#'    using the statistical thresholds.
#' @param adjp_cutoff,p_cutoff,fold_cutoff,mgm_cutoff,ave_cutoff `numeric`
#'    values representing the appropriate statistical threshold,
#'    or `NULL` when a threshold should not be applied.
#' @param int_adjp_cutoff,int_p_cutoff,int_fold_cutoff `numeric`
#'    thresholds to apply only to interaction contrasts.
#' @param confint `logical` passed to `limma::topTable()`, which defines
#'    whether to return confidence intervals for each log2 fold change.
#' @param use_cutoff_colnames `logical` whether to include the
#'    statistical thresholds abbreviated in the `"hit"` colname,
#'    when `define_hits=TRUE`.
#' @param rename_headers `logical` indicating whether to rename
#'    statistical colnames returned by `limma::topTable()` to the
#'    colnames include the contrast name.
#' @param return_fold `logical` whether to return an additional column
#'    with the signed fold change, see `log2fold_to_fold()`.
#' @param merge_df `logical` indicating whether to merge the final
#'    `data.frame` list into one `data.frame`.
#' @param include_ave_expr `logical` indicating whether to retain
#'    the column `"AveExpr"`. This column can be misleading, especially
#'    if the `mgm` (max group mean) threshold is used when determining
#'    statistical hits. This column is mainly useful in reviewing limma
#'    output, since it uses the `"AveExpr"` values to apply its moderated
#'    variance statistic.
#' @param include_group_means `logical` indicating whether to include each
#'    group mean along with the relevant contrast. These values are
#'    helpful, in that they should exactly represent the reported `logFC`
#'    value. Sometimes it is helpful and comforting to see the exact values
#'    used in that calculation.
#' @param collapse_by_gene `logical` indicating whether to apply
#'    `collapse_stats_by_gene` which chooses one "best" exemplar per gene
#'    when there are multiple rows that represent the same gene.
#' @param rename_contrasts `logical` (inactive) which will in future allow
#'    for automated renaming of contrasts.
#' @param sep `character` string used as a delimiter in certain output
#'    colnames.
#' @param int_grep `character` string used to recognize contrasts which
#'    are considered "interaction contrasts". The default pattern recognizes
#'    any contrasts that contain multiple fold changes, recognized by the
#'    presence of more than one hypen `"-"` in the contrast name.
#' @param verbose `logical` indicating whether to print verbose output.
#' 
#' @export
ebayes2dfs <- function
(lmFit3=NULL,
 lmFit1=NULL,
 define_hits=TRUE,
 adjp_cutoff=0.05, 
 p_cutoff=NULL, 
 fold_cutoff=1.5,
 int_adjp_cutoff=adjp_cutoff, 
 int_p_cutoff=p_cutoff,
 int_fold_cutoff=fold_cutoff,
 mgm_cutoff=NULL,
 ave_cutoff=NULL,
 confint=FALSE,
 use_cutoff_colnames=TRUE,
 rename_headers=TRUE,
 return_fold=TRUE, 
 merge_df=FALSE, 
 include_ave_expr=FALSE, 
 include_group_means=TRUE,
 transform_means=c("none", "exp2signed", "10^"),
 collapse_by_gene=FALSE, 
 rename_contrasts=FALSE,
 sep=" ",
 int_grep="[(].+-.+-.+[)]|-.+-",
 trim_colnames=c("t", "B", "F"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert the lmFit3 results of eBayes() into a list of data.frames
   ##
   ## Note the cutoffFold is normal space fold change, but converted to log2fold to compare with limma output
   ## Note cutoffMaxGroupmean is in whatever units are sent to limma... usually log2 intensity or log2 counts
   ##
   ## collapseByGene=TRUE will attempt to produce per-gene results, using pre-defined logic to select the best
   ## exemplar(s) among multiple probes for the same gene.
   ##    1. choose statistical hits first
   ##       a. if multiple hits, same direction, choose them all.
   ##       b. if multiple hits, diff direction,
   ##          i.  mark this gene as multi-direction
   ##          ii. choose best P-value, then take hits with same direction
   ##    2. if no statistical hits, choose entry(ies) above maxMean signal cutoff
   ##    3. If no hits, and no entries above signal cutoff, choose all entries
   ##
   ## renameContrasts=TRUE will rename a fully described contrast into a more human-readable
   ## contrast, e.g.
   ## from: dHSA10GR_EtOH-SW13GR_EtOH
   ## to:   dHSA10GR-SW13GR (EtOH)
   ##
   ## intCutoffPVal and intCutoffAdjPVal are intended to allow
   ## applying a different P-value threshold for interaction effects.
   ##
   ## confint is used by topTable(), when FALSE no confidence intervals are
   ## calculated; when TRUE, by default it calculates 0.95 confidence
   ## intervals, reported as logFC upper and lower bounds, CI.L and CI.R,
   ## respectively.
   ##
   ## Optionally transform the AveExpr values, most useful when reporting normal space values which are stored in log space
   transform_means <- match.arg(transform_means);
   
   ## cutoffMaxGroupMean requires lmFit1
   if (!define_hits) {
      mgm_cutoff <- NULL;
   }
   if ((length(jamba::rmNA(mgm_cutoff)) > 0 || include_group_means) && length(lmFit1) == 0) {
      #stop("To use cutoffMaxGroupMean, lmFit1 must be supplied, from which the group means are obtained.");
      jamba::printDebug("ebayes2dfs(): ",
         c("lmFit1 is required for mgm_cutoff or include_group_means. Setting ",
            "mgm_cutoff=NULL",
            ", and ",
            "include_group_means=FALSE"),
         sep="");
      include_group_means <- FALSE;
      mgm_cutoff <- NULL;
   }

   ## TODO: Allow some curation of labels here
   contrastNames <- colnames(coef(lmFit3));
   contrastLabels <- jamba::nameVector(contrastNames, contrastNames);
   is_interaction <- grepl(int_grep,
      ignore.case=TRUE,
      contrastNames);
   
   if (define_hits) {
      ## Add cutoff parameters to the colnames, optional
      cutoff_l <- list(
         mgm=mgm_cutoff,
         p=p_cutoff,
         adjp=adjp_cutoff,
         fc=fold_cutoff,
         intp=int_p_cutoff,
         intadjp=int_adjp_cutoff,
         intfc=int_fold_cutoff);
      cutoff_df <- as.data.frame(jamba::rmNULL(cutoff_l));
      if (!any(is_interaction) && jamba::igrepHas("int", colnames(cutoff_df))) {
         cutoff_df <- cutoff_df[,jamba::unvigrep("int", colnames(cutoff_df)),drop=FALSE];
      }
      if (nrow(cutoff_df) == 0) {
         define_hits <- FALSE;
      }
   }
   if (define_hits) {
      cutoff_string_df <- cutoff_df;
      for (i in colnames(cutoff_df)) {
         cutoff_string_df[[i]] <- paste0(i, cutoff_df[[i]]);
         ## if column is "int"
         ## remove if values are all identical to "non-int" cutoff
         if (jamba::igrepHas("^int", i)) {
            j <- gsub("^int", "", i);
            if (j %in% colnames(cutoff_df)) {
               if (all(cutoff_df[[i]] == cutoff_df[[j]])) {
                  cutoff_string_df[,i] <- list(NULL);
               }
            }
         }
      }
      if (any(is_interaction)) {
         int_use <- jamba::provigrep(
            c("mgm", "ave", "^intadjp", "^intp", "^intfc", "^int"),
            colnames(cutoff_string_df))
         if (length(int_use) > 0) {
            int_cutoff_string <- jamba::pasteByRow(
               cutoff_string_df[,int_use,drop=FALSE],
               sep=sep);
         } else {
            int_cutoff_string <- rep("", nrow(cutoff_df));
         }
      }
      nonint_use <- jamba::unvigrep("^int",
         jamba::provigrep(
            c("mgm", "ave", "^adjp", "^p", "^fc"),
            colnames(cutoff_string_df)));
      if (length(nonint_use) > 0) {
         cutoff_string <- jamba::pasteByRow(
            cutoff_string_df[,nonint_use,drop=FALSE],
            sep=sep);
      } else {
         cutoff_string <- rep("", nrow(cutoff_df));
      }
      if (verbose) {
         jamba::printDebug("ebayes2dfs(): ",
            "cutoff_df:");
         print(cutoff_df);
         jamba::printDebug("ebayes2dfs(): ",
            "cutoff_string_df:");
         print(cutoff_string_df);
         jamba::printDebug("ebayes2dfs(): ",
            "cutoff_string:",
            cutoff_string);
      }
      #return(cutoff_string_df);
   }

   ## assign rownames if not present in lmFit1$coefficients
   if (length(lmFit1) > 0 && length(rownames(lmFit1$coefficients)) == 0) {
      jamba::printDebug("ebayes2dfs(): ",
         c("Note there are no",
            " rownames(lmFit1$coefficients)"),
         sep="",
         fgText=c("darkorange1", "red1"));
      if ("genes" %in% names(lmFit1)) {
         rownames(lmFit1$coefficients) <- rownames(lmFit1$genes);
      }
   }
   if (length(rownames(lmFit3$coefficients)) == 0) {
      jamba::printDebug("ebayes2dfs(): ",
         c("Note there are no",
            " rownames(lmFit3$coefficients)"),
         sep="",
         fgText=c("darkorange1", "red1"));
      if ("genes" %in% names(lmFit3)) {
         rownames(lmFit3$coefficients) <- rownames(lmFit3$genes);
      }
   }

   ## TODO: handle cases with zero residual degrees of freedom,
   ## where we would not have a P-value but still have fold changes.
   ## Examples would be per-patient fold changes.
   ##
   ## lmTopTables is a list:
   ## - named by contrastNames
   ## - containing elements "iTopTableDF"
   ## - if collapseByGene=TRUE
   ##    - "iTopTableByGene"
   ##    - "multiDirProbes" 
   lmTopTables <- lapply(jamba::nameVector(contrastNames), function(i){
      retVals <- list();
      iLabel <- contrastLabels[i];
      ## Detect whether the contrast is an interaction effect (2-way or 3-way)
      ## or is a simple pairwise style t-test
      isInteraction <- jamba::igrepHas(int_grep, iLabel);
      if (verbose) {
         if (isInteraction) {
            jamba::printDebug("ebayes2dfs(): ",
               c("Interaction effect detected for contrast:",
                  iLabel),
               sep="",
               fgText=c("darkorange1", "purple"));
         } else {
            jamba::printDebug("ebayes2dfs(): ",
               c("Evaluting contrast:",
                  iLabel),
               sep="");
         }
      }
      #if (rename_contrasts) {
      #   iLabel <- renameTwoFactorComparisons(iLabel);
      #}
      
      ## Assemble top table, handling single replicate data in a specific way
      if (!any(lmFit3$df.residual > 0)) {
         jamba::printDebug("ebayes2dfs(): ",
            "No values for df.residual>0.",
            fgText=c("darkorange1", "red1"));
         if (confint) {
            iTopTable <- data.frame(
               check.names=FALSE,
               stringsAsFactors=FALSE,
               lmFit3$genes,
               logFC=lmFit3$coefficients[,i],
               CI.L=lmFit3$coefficients[,i],
               CI.R=lmFit3$coefficients[,i],
               adj.P.Val=1,
               P.Value=1,
               AveExpr=lmFit3$Amean);
         } else {
            iTopTable <- data.frame(
               check.names=FALSE,
               stringsAsFactors=FALSE,
               lmFit3$genes,
               logFC=lmFit3$coefficients[,i],
               adj.P.Val=1,
               P.Value=1,
               AveExpr=lmFit3$Amean);
         }
      } else {
         iTopTable <- limma::topTable(
            lmFit3,
            coef=i,
            sort.by="none",
            number=nrow(lmFit3),
            confint=confint);
      }
      ## Optionally remove extraneous colnames
      if (length(trim_colnames) > 0) {
         iTopTable <- iTopTable[,setdiff(colnames(iTopTable), trim_colnames), drop=FALSE];
      }
      
      ## Optionally include group mean and maxgroupmean values
      if (length(lmFit1) > 0) {
         ## Improve maxGroupMean by using the specific groups included in the contrast
         iCoefCols <- names(which(lmFit3$contrasts[,i] != 0));
         iCoefLabs <- paste(iCoefCols,
            sep=sep,
            "mean");
         coef_match <- match(rownames(iTopTable),
            rownames(lmFit1$coefficients));
         gm_m <- jamba::renameColumn(
            lmFit1$coefficients[coef_match, iCoefCols, drop=FALSE],
            from=iCoefCols,
            to=iCoefLabs);
         mgm <- matrixStats::rowMaxs(gm_m,
            na.rm=TRUE);
         if (include_group_means) {
            iTopTable[,colnames(gm_m)] <- gm_m;
         }
         iTopTable[,"mgm"] <- mgm;
      }
      
      ## Apply statistical thresholds
      if (define_hits || collapse_by_gene) {
         probe_colname <- head(colnames(iTopTable), 1);
         if (collapse_by_gene) {
            if (verbose) {
               jamba::printDebug("ebayes2dfs(): ",
                  "collapse_by_gene.");
            }
            gene_colname <- head(jamba::provigrep(
               c("GeneName", "geneSymbol", "gene"),
               colnames(iTopTable)), 1);
            isGenes <- nameVector(iTopTable[,gene_colname],
               rownames(iTopTable));
         }
         if (isInteraction) {
            pcol <- "intp";
            adjpcol <- "intadjp";
            foldcol <- "intfc";
         } else {
            pcol <- "p";
            adjpcol <- "adjp";
            foldcol <- "fc";
         }
         
         ## iterate each set of thresholds by row in cutoff_df
         ## but only iterate unique set of cutoff thresholds
         if (isInteraction) {
            k_rows <- match(unique(int_cutoff_string),
               int_cutoff_string);
         } else {
            k_rows <- match(unique(cutoff_string),
               cutoff_string);
         }
         for (k in k_rows) {
            if (isInteraction) {
               hit_colname <- paste("hit",
                  sep=sep,
                  int_cutoff_string[k]);
            } else {
               hit_colname <- paste("hit",
                  sep=sep,
                  cutoff_string[k]);
            }
            if (verbose) {
               jamba::printDebug("ebayes2dfs(): ",
                  c("hit_colname:",
                     hit_colname),
                  sep="");
            }
            ## Call utility function to get hit flags
            mgm_cutoff <- ifelse("mgm" %in% colnames(cutoff_df),
               cutoff_df[k,"mgm"],
               NA);
            ave_cutoff <- ifelse("ave" %in% colnames(cutoff_df),
               cutoff_df[k,"ave"],
               NA);
            p_cutoff <- ifelse(pcol %in% colnames(cutoff_df),
               cutoff_df[k,pcol],
               NA);
            adjp_cutoff <- ifelse(adjpcol %in% colnames(cutoff_df),
               cutoff_df[k,adjpcol],
               NA);
            fold_cutoff <- ifelse(foldcol %in% colnames(cutoff_df),
               cutoff_df[k,foldcol],
               NA);
            hit_values <- mark_stat_hits(x=iTopTable,
               adjp_cutoff=adjp_cutoff,
               p_cutoff=p_cutoff,
               fold_cutoff=fold_cutoff,
               mgm_cutoff=mgm_cutoff,
               ave_cutoff=ave_cutoff,
               adjp_colname="adj.P.Val",
               p_colname="P.Value",
               logfc_colname="logFC",
               mgm_colname="mgm",
               ave_colname="AveExpr");
            iTopTable[[hit_colname]] <- hit_values;
         }
      }
      
      ## Optionally transform intensities, e.g. exponentiating log2 values
      ## Note that the 2^ conversion now subtracts 1, since the typical
      ## transformation is log2(1+x)
      if (1 == 2) {
         if (transformAveExpr %in% c("2^")) {
            #iTopTable[,"AveExpr"] <- 2^iTopTable[,"AveExpr"];
            iTopTable[,"AveExpr"] <- 2^iTopTable[,"AveExpr"] - 1;
            if ("maxGroupMean" %in% colnames(iTopTable)) {
               #iTopTable[,"maxGroupMean"] <- 2^iTopTable[,"maxGroupMean"];
               iTopTable[,"maxGroupMean"] <- 2^iTopTable[,"maxGroupMean"] - 1;
            }
         } else if (transformAveExpr %in% c("10^")) {
            #iTopTable[,"AveExpr"] <- 10^iTopTable[,"AveExpr"];
            iTopTable[,"AveExpr"] <- 10^iTopTable[,"AveExpr"] - 1;
            if ("maxGroupMean" %in% colnames(iTopTable)) {
               #iTopTable[,"maxGroupMean"] <- 10^iTopTable[,"maxGroupMean"];
               iTopTable[,"maxGroupMean"] <- 10^iTopTable[,"maxGroupMean"] - 1;
            }
         }
      }
      ## Optionally convert log2 fold change to normal fold change
      if (return_fold) {
         iTopTable[,"fold"] <- log2fold_to_fold(iTopTable[,"logFC"]);
      }
      ################################################
      ## Per-gene collapse:
      ##
      ## Note we ultimately define per-gene rows using the first stats hit criteria,
      ## otherwise too many columns will be created.  E.g. changing the stats hit criteria
      ## changes the prioritization of probes to include in the per-gene row, so it affects
      ## the fold change, groupMean, P-Value, and adj-P-Val.
      ## The current solution is to use only the first hit filters, so only one set of these
      ## values is propagated downstream.
      ## However, the "hit" columns can represent multiple stats hit criteria.
      ##
      ## TODO: implement the interaction P-value cutoffs here as well
      if (collapse_by_gene) {
         if (verbose) {
            printDebug("collapseByGene iTopTable:");
            print(head(iTopTable));
            printDebug("geneColname:", geneColname);
         }
         iTopTableByGeneL <- collapseTopTableByGene(iTopTable,
            geneColname=geneColname,
            aveExprColname="AveExpr",
            maxGroupMeanColname="maxGroupMean",
            adjPvalColname="adj.P.Val",
            PvalueColname="P.Value",
            logFCcolname="logFC",
            cutoffAveExpr=cutoffAveExpr[1],
            cutoffMaxGroupMean=cutoffMaxGroupMean[1],
            cutoffAdjPVal=cutoffAdjPVal[1],
            cutoffPVal=cutoffPVal[1],
            cutoffFold=cutoffFold[1]);
         iTopTableByGene <- iTopTableByGeneL$iTopTableByGene;
         retVals$top_bygene_df <- iTopTableByGene;
         retVals$multidir_probes <- iTopTableByGeneL$multiDirProbes;
      }
      ## Optionally rename column headers to include the contrast name
      gene_colnames <- intersect(colnames(iTopTable),
         c(colnames(lmFit3$genes), probe_colname));
      if (rename_headers) {
         ## Note we do rename maxGroupMean now,
         ## since we calculate it per each specific contrast
         rename_from <- setdiff(
            jamba::unvigrep("AveExpr| mean$|gene|probe",
               colnames(iTopTable)),
            gene_colnames);
         rename_to <- paste(rename_from,
            iLabel,
            sep=sep);
         iTopTable <- jamba::renameColumn(iTopTable,
            from=rename_from,
            to=rename_to);
      }
      if (!include_ave_expr) {
         iTopTable <- iTopTable[,jamba::unvigrep("AveExpr", colnames(iTopTable)),drop=FALSE];
      }
      rownames(iTopTable) <- jamba::makeNames(iTopTable[,1]);
      retVals$top_df <- iTopTable;
      retVals;
   });
   
   ## Now prepare the data to return
   if (merge_df) {
      lmTopTablesAll <- jamba::mergeAllXY(lapply(lmTopTables, function(i){
         i$top_df;
      }));
      if (collapse_by_gene) {
         lmTopTablesAllG <- jamba::mergeAllXY(lapply(lmTopTables, function(i){
            i$top_bygene_df;
         }));
      }
      ## Re-order columns so the genes, then hits, appear first
      colname_order <- jamba::provigrep(c(gene_colnames, "^hit", "."),
         colnames(lmTopTablesAll));
      lmTopTablesAll <- lmTopTablesAll[,colname_order, drop=FALSE];
      rownames(lmTopTablesAll) <- jamba::makeNames(lmTopTablesAll[,head(gene_colnames, 1)]);
   } else {
      lmTopTablesAll <- lapply(lmTopTables, function(i){
         i$top_df;
      });
      if (collapse_by_gene) {
         lmTopTablesAllG <- lapply(lmTopTables, function(i){
            i$top_bygene_df;
         });
      }
   }
   
   if (define_hits) {
      attr(lmTopTablesAll, "cutoff_df") <- cutoff_df;
   }
   if (collapse_by_gene) {
      if (define_hits) {
         attr(lmTopTablesAllG, "cutoff_df") <- cutoff_df;
      }
      retList <- list(top_df=lmTopTablesAll,
         top_bygene_df=lmTopTablesAllG);
      return(retList);
   }
   return(lmTopTablesAll);
}

