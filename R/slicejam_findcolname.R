#' Find colname by string or pattern
#' 
#' Find colname by string or pattern
#' 
#' This function is a simple utility function intended to help
#' find the most appropriate matching colname given one or more
#' character strings or patterns.
#' 
#' It returns the first best matching result, but can return
#' multiple results in order of preference if `max=Inf`.
#' 
#' The order of matching:
#' 
#' 1. Match the exact colname.
#' 2. Match case-insensitive colname.
#' 3. Match the beginning of each colname.
#' 4. Match the end of each colname.
#' 5. Match anywhere in each colname.
#' 
#' The goal is to use something like `c("p.value", "pvalue", "pval")`
#' and be able to find colnames with these variations:
#' 
#' * `P.Value`
#' * `P.Value Group-Control`
#' * `Group-Control P.Value`
#' * `pvalue`
#' 
#' Even if the data contains `c("P.Value", "adj.P.Val")` as returned
#' by `limma::topTable()` for example, the pattern `c("p.val")` will
#' preferentially match `"P.Value"` and not `"adj.P.Val"`.
#' 
#' @param pattern `character` vector of text strings and/or regular
#'    expression patterns.
#' @param x `data.frame` or other object that contains `colnames(x)`.
#' @param max `integer` maximum number of entries to return.
#' @param index `logical` indicating whether to return the column index,
#'    that is the column number.
#' @param require_non_na `logical` indicating whether to require at
#'    least one non-`NA` value in the matching colname. When
#'    `require_non_na=TRUE` and all values in a column are `NA`,
#'    that colname is not returned by this function.
#' @param exclude_pattern `character` vector of colnames or patterns
#'    to exclude from returned results.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' x <- data.frame(
#'    `Gene`=paste0("gene", LETTERS[1:25]), 
#'    `log2fold Group-Control`=rnorm(25)*2,
#'    `P.Value Group-Control`=10^-rnorm(25)^2,
#'    check.names=FALSE);
#' x[["fold Group-Control"]] <- log2fold_to_fold(x[["log2fold Group-Control"]]);
#' x[["adj.P.Val Group-Control"]] <- x[["P.Value Group-Control"]];
#' 
#' print(head(x));
#' find_colname(c("p.val", "pval"), x);
#' find_colname(c("fold", "fc", "ratio"), x);
#' find_colname(c("logfold", "log2fold", "lfc", "log2ratio", "logratio"), x);
#' 
#' ## use exclude_pattern
#' ## if the input data has no "P.Value" but has "adj.P.Val"
#' y <- x[,c(1,2,4,5)];
#' print(head(y));
#' find_colname(c("p.val"), y, exclude_pattern=c("adj"))
#' 
find_colname <- function
(pattern,
 x,
 max=1,
 index=FALSE,
 require_non_na=TRUE,
 col_types=NULL,
 exclude_pattern=NULL,
 verbose=FALSE,
 ...)
{
   ##
   if (length(pattern) == 0) {
      return(NULL);
   }
   x_colnames <- colnames(x);
   if (length(x_colnames) == 0) {
      return(NULL);
   }

   ## col_types
   if (length(col_types) > 0 && jamba::igrepHas("data.frame|tbl|data.table", class(x))) {
      col_classes <- sapply(colnames(x), function(i){
         class(x[[i]])
      });
      x_keep <- (col_classes %in% col_types);
      if (verbose && any(!x_keep)) {
         jamba::printDebug("find_colname(): ",
            "applied col_types filter and removed:",
            x_colnames[!x_keep]);
      }
      x_colnames <- x_colnames[x_keep];
   }
   
   ## require_non_na
   if (require_non_na) {
      x_colnames <- x_colnames[sapply(x_colnames, function(icol){
         any(!is.na(x[[icol]]))
      })]
   }
   ## if no colnames remain, return NULL
   if (length(x_colnames) == 0) {
      return(x_colnames);
   }

   ## Optional exclude_pattern
   if (length(exclude_pattern) > 0) {
      exclude_colnames <- find_colname(pattern=exclude_pattern,
         x=x,
         max=Inf,
         index=FALSE,
         require_non_na=FALSE,
         exclude_pattern=NULL,
         verbose=FALSE);
      if (length(exclude_colnames) > 0) {
         if (verbose) {
            jamba::printDebug("find_colname(): ",
               "exclude_colnames:",
               exclude_colnames);
         }
         x_colnames <- setdiff(x_colnames,
            exclude_colnames);
      }
   } else {
      exclude_colnames <- NULL;
   }
   
   start_pattern <- paste0("^", pattern);
   end_pattern <- paste0(pattern, "$");
   
   if (any(pattern %in% x_colnames)) {
      ## 1. max exact colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning exact match.");
      }
      x_vals <- intersect(pattern, x_colnames);
   } else if (any(tolower(pattern) %in% tolower(x_colnames))) {
      ## 2. max exact colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning exact case-insensitive match.");
      }
      x_match <- jamba::rmNA(match(tolower(pattern), tolower(x_colnames)));
      x_vals <- x_colnames[x_match];
   } else if (jamba::igrepHas(paste(collapse="|", start_pattern), x_colnames)) {
      ## 3. match start of each colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning match to colname start.");
      }
      x_vals <- unique(jamba::provigrep(start_pattern, x_colnames));
   } else if (jamba::igrepHas(paste(collapse="|", end_pattern), x_colnames)) {
      ## 4. match end of each colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning match to colname end.");
      }
      x_vals <- unique(jamba::provigrep(end_pattern, x_colnames));
   } else if (jamba::igrepHas(paste(collapse="|", pattern), x_colnames)) {
      ## 5. match any part of each colname
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "Returning match to part of colname.");
      }
      x_vals <- unique(jamba::provigrep(pattern, x_colnames));
      return(head(x_vals, max));
   } else {
      if (verbose) {
         jamba::printDebug("find_colname(): ",
            "No match found.");
      }
      x_vals <- NULL;
   }
   
   if (index && length(x_vals) > 0) {
      x_vals <- unique(match(x_vals, x_colnames));
   }
   return(head(x_vals, max));
}
