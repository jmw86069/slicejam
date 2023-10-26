
#' Sort contrasts by factor
#' 
#' Sort contrasts by factor
#' 
#' @param x `character` vector of contrast names
#' @param contrast_style `character` string to define the return type:
#'    * `"contrast"` returns full contrast names
#'    * `"comp"` returns abbraviate comp as from `jamses::contrast2comp()`
#' @param ... additional arguments are passed to
#'    `jamses::contrasts_to_factors()`
#' 
#' @export
sort_contrast_names <- function
(x,
 factor_names=NULL,
 contrast_style=c("contrast",
    "comp"),
 verbose=FALSE,
 ...)
{
   # validate arguments
   contrast_style <- match.arg(contrast_style);
   #
   factors_df <- jamses::contrasts_to_factors(x)
   if (length(factor_names) == ncol(factors_df)) {
      colnames(factors_df) <- factor_names
   }
   
   factorcols1 <- colnames(factors_df)
   factorcols <- c(tail(factorcols1, -1), head(factorcols1, 1))
   # convert each column to factor so we maintain any implicit ordering
   for (icol in factorcols) {
      factors_df[[icol]] <- factor(factors_df[[icol]],
         levels=unique(factors_df[[icol]]));
   }
   # define depth of the comparison (how many factors have a comparison)
   factors_df$depth <- rowSums(do.call(cbind,
      lapply(factorcols, function(icol){
         grepl("-", factors_df[[icol]]) * 1
      })
   ))
   # store which factor was compared
   factors_df$compared_factors <- jamba::pasteByRow(do.call(cbind,
      lapply(factorcols, function(icol){
         ifelse(grepl("-", factors_df[[icol]]), icol, "")
      })),
      sep=":")
   # convert to factor to main ordering
   factors_df$compared_factors <- factor(factors_df$compared_factors,
      levels=provigrep(c(paste0("^", factorcols, "$"), factorcols, "."),
         unique(factors_df$compared_factors)))
   # store the comparisons used, so they can inform the sort order
   factors_df$sort_field <- jamba::pasteByRow(do.call(cbind,
      lapply(factorcols, function(icol){
         ifelse(grepl("-", factors_df[[icol]]), as.character(factors_df[[icol]]), "")
      })),
      sep=":")
   
   # apply the sort
   sorted_factors_df <- jamba::mixedSortDF(factors_df,
      byCols=c("depth",
         "compared_factors",
         "sort_field",
         factorcols))
   if (verbose) {
      jamba::printDebug("sort_contrast_names(): ",
         "sorted_factors_df:");
      print(sorted_factors_df);
   }
   sorted_comps <- jamba::pasteByRow(
      sorted_factors_df[, factorcols1, drop=FALSE],
      sep=":");
   if ("comp" %in% contrast_style) {
      return(unname(sorted_comps))
   }
   sorted_contrast_names <- jamses::comp2contrast(sorted_comps);
   names(sorted_contrast_names) <- sorted_comps;
   return(sorted_contrast_names);
}
