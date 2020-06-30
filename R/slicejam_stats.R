
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
         if (length(genes) < nrow(se) || length(samples) < ncol(se)) {
            assays(se)[[output_assay_name]] <- assays(se)[[assay_name]];
            assays(se)[[output_assay_name]][] <- NA;
         }
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
