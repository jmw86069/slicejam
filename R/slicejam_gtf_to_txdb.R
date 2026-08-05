
#' Convert GTF to TxDb with caching
#'
#' Convert GTF to TxDb with caching, called by other
#' SliceJam functions.
#' 
#' * It accepts GTF as a filename, or web URL.
#' * Default `save_txdb=TRUE` will try to save to the current
#' working directory, using argument `txdb_path`.
#' Or you can specify a specific filename `txdb_file`.
#' * When `txdb_file` is supplied without any slashes,
#' it uses `txdb_path`. If `txdb_file` contains a slash,
#' then it is used without `txdb_path`.
#' 
#' @param gtf `character` path or URL to a GTF or GFF3 file;
#'    or `GRanges` object suitable for use with
#'    `txdbmaker::makeTxDbFromGRanges()` to create a `TxDb`.
#' @param save_txdb `logical` default TRUE, whether to save the
#'    `TxDb` data to a `'.txdb'` file for future re-use.
#' @param force_refresh `logical` default FALSE, whether to
#'    force the `TxDb` to be created even if there exists a
#'    corresponding `'.txdb'` file.
#' @param txdb_file `character` name of a file. When the file
#'    does not contain any slash `'/'` character, it will be
#'    used together with `txdb_path` to create the file path.
#' @param txdb_path `character` default `'.'` used to form the
#'    full path to the output `'.txdb'` file. When `txdb_file`
#'    contains a full path, containing any `'/'` slash character,
#'    this argument is ignored.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are passed to the corresponding
#'    `txdbmaker::makeTxDb()` function.
#' 
#' @returns `TxDb` object as defined in Bioconductor package
#' `txdbmaker::makeTxDb`.
#' 
#' @keywords internal
#' @noRd
convert_gtf_to_txdb <- function
(gtf=NULL,
 save_txdb=TRUE,
 force_refresh=FALSE,
 txdb_file=NULL,
 txdb_path=".",
 verbose=FALSE,
 ...)
{
   #
   if (length(gtf) == 0) {
      stop("'gtf' is required.")
   }

   define_txdb_filepath <- function(txdb_file, txdb_path) {
      if (length(txdb_file) == 0 || nchar(txdb_file) == 0) {
         return(NULL)
      }
      if (grepl("[\\/]", txdb_file)) {
         return(txdb_file);
      }
      if (length(txdb_path) == 0 || nchar(txdb_path) == 0) {
         txdb_path <- "."
      }
      file.path(txdb_file, txdb_path)
      #
   }
   txdb_filepath <- NULL;

   ###########################################
   ## Accept GRanges
   if (inherits(gtf, "GRanges")) {
      if (isTRUE(save_txdb) && length(txdb_file) == 0) {
         stop("'txdb_file' is required when supplying GRanges.")
      } else {
         txdb_filepath <- define_txdb_filepath(txdb_file, txdb_path)
      }
      if (
         isTRUE(save_txdb) &&
            file.exists(txdb_filepath) &&
            !isTRUE(force_refresh)
      ) {
         # loadDb from TxDb file
         if (verbose) {
            jamba::printDebug(
               "convert_gtf_to_txdb(): ",
               "Loading existing txdb file:",
               txdb_filepath
            )
         }
         refgene_txdb <- AnnotationDbi::loadDb(txdb_filepath)
         #
      } else {
         if (verbose) {
            jamba::printDebug(
               "convert_gtf_to_txdb(): ",
               "Creating TxDb from GRanges."
            )
         }
         refgene_txdb <- jamba::call_fn_ellipsis(
            txdbmaker::makeTxDbFromGRanges,
            gr = gtf,
            ...
         )
         if (isTRUE(save_txdb)) {
            # save TxDb to file
            AnnotationDbi::saveDb(refgene_txdb, file = txdb_filepath)
            if (verbose) {
               jamba::printDebug(
                  "convert_gtf_to_txdb(): ",
                  "Saved txdb file:",
                  txdb_filepath
               )
            }
            #
         }
      }
      #
   } else {
      ###########################################
      ## GTF file or URL
      gtf <- head(gtf, 1);
      if (!(
            grepl("^(ftp|http|sftp|ssh)://", gtf) ||
            file.exists(gtf))) {
         stop(paste0("'gtf' is not a URL and file not found:", gtf))
      }
      # Todo: Accept URL, and clean URL with parameters: file.gtf?format=3
      # remove the extension: .gtf, .gff, .gff3,
      # .gtf.gz, .gtf.xz, .gtf.Z, .gtf.zip
      gtf_base <- gsub(
         "[.](gff|gtf|gff3)([-_.](zip|gz|xz|Z)|)$",
         "",
         ignore.case=TRUE,
         gsub("[?].*$", "",
            basename(gtf)))

      # txdb_file stores previous conversion GTF->txdb
      txdb_file <- paste0(gtf_base, ".txdb")
      txdb_filepath <- define_txdb_filepath(txdb_file, txdb_path);
      refgene_txdb <- NULL
      # If txdb_file exists, load it
      if (file.exists(txdb_filepath)) {
         if (verbose) {
            jamba::printDebug(
               "convert_gtf_to_txdb(): ",
               "Loading existing txdb file:",
               txdb_filepath
            )
         }
         refgene_txdb <- AnnotationDbi::loadDb(txdb_filepath)
      } else {
         #
         ## 15-20 seconds for human GTF
         if (verbose) {
            jamba::printDebug(
               "convert_gtf_to_txdb(): ",
               "Creating txdb from gtf:",
               gtf
            )
         }
         refgene_txdb <- jamba::call_fn_ellipsis(
            txdbmaker::makeTxDbFromGFF,
            file = gtf,
            ...
         )
         if (isTRUE(save_txdb)) {
            # save TxDb to file
            AnnotationDbi::saveDb(refgene_txdb,
               file = txdb_filepath)
            if (verbose) {
               jamba::printDebug(
                  "convert_gtf_to_txdb(): ",
                  "Saved txdb file:",
                  txdb_filepath
               )
            }
            #
         }
      }
   }
   return(refgene_txdb)
}

#' Convert GTF to tx2gene with caching
#'
#' Convert GTF to tx2gene with caching, called by other
#' SliceJam functions. It relies on `splicejam::makeTx2geneFromGtf()`
#' and as such, the '...' arguments are passed through.
#'
#' @returns `data.frame` with columns: 'transcript_id', 'gene_name'
#' @keywords internal
#' @noRd
convert_gtf_to_tx2gene <- function(
   gtf = NULL,
   save_tx2gene = TRUE,
   force_refresh = FALSE,
   tx2gene_file = NULL,
   tx2gene_path = ".",
   verbose = FALSE,
   ...
) {
   #
   if (length(gtf) == 0) {
      stop("'gtf' is required.")
   }
   gtf <- head(gtf, 1)
   if (
      !(grepl("^(ftp|http|sftp|ssh)://", gtf) ||
         file.exists(gtf))
   ) {
      stop(paste0("'gtf' is not a URL and file not found:", gtf))
   }

   define_tx2gene_filepath <- function(tx2gene_file, tx2gene_path) {
      if (length(tx2gene_file) == 0 || nchar(tx2gene_file) == 0) {
         return(NULL)
      }
      if (grepl("[\\/]", tx2gene_file)) {
         return(tx2gene_file)
      }
      if (length(tx2gene_path) == 0 || nchar(tx2gene_path) == 0) {
         tx2gene_path <- "."
      }
      file.path(tx2gene_file, tx2gene_path)
      #
   }
   tx2gene_filepath <- NULL

   # Todo: Accept URL, and clean URL with parameters: file.gtf?format=3
   # remove the extension: .gtf, .gff, .gff3,
   # .gtf.gz, .gtf.xz, .gtf.Z, .gtf.zip
   gtf_base <- gsub(
      "[.](gff|gtf|gff3)([-_.](zip|gz|xz|Z)|)$",
      "",
      ignore.case = TRUE,
      gsub("[?].*$", "", basename(gtf))
   )

   # txdb_file stores previous conversion GTF->txdb
   tx2gene_file <- paste0(gtf_base, ".tx2gene.txt")
   tx2gene_filepath <- define_tx2gene_filepath(tx2gene_file, tx2gene_path)
   tx2geneDF <- NULL
   # If txdb_file exists, load it
   if (
      !isTRUE(force_refresh) &&
         isTRUE(save_tx2gene) &&
         file.exists(tx2gene_filepath)
   ) {
      if (verbose) {
         jamba::printDebug(
            "convert_gtf_to_tx2gene(): ",
            "Loading existing tx2gene file:",
            tx2gene_filepath
         )
      }
      tx2geneDF <- data.table::fread(
         txdb_filepath,
         sep = "\t",
         data.table = FALSE
      )
      # Todo: validate colnames
   } else {
      if (verbose) {
         jamba::printDebug(
            "convert_gtf_to_tx2gene(): ",
            "Creating tx2gene from GTF:",
            gtf
         )
      }
      tx2geneDF <- splicejam::makeTx2geneFromGtf(
         gtf,
         geneAttrNames = geneAttrNames,
         txAttrNames = txAttrNames,
         geneFeatureType = geneFeatureType,
         verbose = verbose,
         txFeatureType = txFeatureType,
         ...
      )
      #
      if (isTRUE(save_tx2gene)) {
         data.table::fwrite(
            tx2geneDF,
            file = tx2gene_file,
            quote = FALSE,
            sep = "\t"
         )
         if (verbose) {
            jamba::printDebug(
               "convert_gtf_to_tx2gene(): ",
               "Saved tx2gene file:",
               tx2gene_filepath
            )
         }
         #
      }
   }


   #
   if (!file.exists(tx2gene_file) && file.exists(basename(tx2gene_file))) {
      tx2gene_file <- basename(tx2gene_file)
   }
   if (file.exists(tx2gene_file)) {
      tx2geneDF <- data.table::fread(
         tx2gene_file,
         sep = "\t",
         data.table = FALSE
      )
   } else {
      tx2geneDF <- splicejam::makeTx2geneFromGtf(
         gtf,
         geneAttrNames = geneAttrNames,
         txAttrNames = txAttrNames,
         geneFeatureType = geneFeatureType,
         verbose = verbose,
         txFeatureType = txFeatureType,
         ...
      )
      # save to a file
      tryCatch(
         {
            data.table::fwrite(tx2geneDF, file = tx2gene_file, sep = "\t")
         },
         error = function(e) {
            data.table::fwrite(
               tx2geneDF,
               file = basename(tx2gene_file),
               sep = "\t"
            )
         }
      )
   }
   # Todo: validate colnames?

   return(tx2geneDF)
}

