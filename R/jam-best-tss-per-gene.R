
#' Define the best TSS per gene from 'Salmon' quant data
#'
#' Define the best TSS per gene from 'Salmon' quant data,
#' using transcript quantitation produced by `tximport::tximport()`
#' to determine the TSS and TTS with the highest supporting
#' signal.
#'
#' 'Salmon quant' produces estimated quantitation for
#' individual transcript isoforms, whose data are most often
#' imported using `tximport::tximport()`. As such, similar
#' tools such as 'kallisto' are also capable of providing
#' equivalent input data.
#'
#' Abundance estimates include transcripts in the corresponding
#' transcriptome index file used by 'Salmon' or 'kallisto'.
#' For the purpose of the workflow used here, the data are
#' imported using `tximport::tximport()`.
#' We typically use `countsFromAbundance='lengthScaledTPM'`
#' however, the abundance itself can be used and not the
#' estimated pseudo-counts.
#'
#' The workflow essentially re-uses `tximport::summarizeByGene()`
#' where 'transcript_id' values are grouped by gene and transcript
#' start site (TSS) to aggregate abundance estimates by TSS.
#' Then the highest abundance TSS is selected as "preferred".
#' 
#' # Handling both TSS and TTS
#' 
#' When providing only the 'gtf', it will calculate all TSS
#' and TTS sites, and both will be processed.
#' 
#' When both TSS and TTS sites are provided, it performs the
#' analysis in tandem:
#' 
#' 1. It determines the 'best TSS' site for each gene.
#' 2. It identifies all transcripts which have these TSS sites.
#' 3. It subsets TTS 'tx2tts' to use only these transcripts,
#' so that the 'best TTS' site will only be evaluated using
#' transcripts compatible with the best 'TSS' site.
#' 
#' There are rare instances where the 'best TTS' site may
#' be upstream from the 'best TSS' site, which does not make
#' sense. While this may be true, that a particular TTS has
#' more signal than a downstream TSS site annotated to the
#' same gene locus, for the purpose of this method, the
#' TTS is assumed to represent transcription that originated
#' from the 'best TSS' site.
#' 
#' To retain the 'best TTS' independent of the 'best TSS',
#' process the TSS `'tx2tss'` and the TTS `'tx2tts'` data
#' in two independent steps.
#'
#' ## Additional Comments
#'
#' * Note that `"transcript_id"` values are not returned.
#' There are often many transcript `'transcript_id'`
#' with the same TSS.
#' The transcripts are not essential to this processing stage.
#' * It is possible to recover the transcripts by reviewing
#' the `'tx2tss'` column 'tss', and comparing with the
#' 'best TSS' for a particular gene.
#' 
#'
#' @family slicejam genome regions
#'
#' @param txiTx `list` or `SummarizedExperiment` with all assay names:
#'    `c("counts", "abundance", "length")`
#' @param groups `character` (default NULL) used to define sample groups,
#'    useful to use the max group mean instead of global mean.
#'    * Using proper `groups` is recommended.
#'    * It can use a column in `colData(txiTx)` when supplying a
#'    `SummarizedExperiment` object.
#' @param tx2tss `data.frame` that must contain two columns:
#'    * `"transcript_id"` - which matches rownames of txiTx.
#'    * `"tss"` - in the exact form `"[gene_name]_[seqnames]_[start]_[strand]"`
#'    for example `"5_8S_rRNA_chr14_16057472"` is recognized to have these
#'    fields: `"5_8S_rRNA", "chr14", "16057472", "+"`
#' @param tx2tts `data.frame` with TTS sites, with the same
#'    format as 'tx2tss'.
#' @param gtf `character` path or URL to a GTF file, used when both
#'    'tx2tss' and 'tx2tts' are NULL.
#'    * In this case, the GTF is used
#'    to define the TSS and TTS sites for each transcript, for each
#'    gene locus.
#'    * The GTF is processed internally by `convert_gtf_to_txdb()`
#'    and `convert_gtf_to_tx2gene()` both of which convert the GTF
#'    to `TxDb`, and either load a previous `'.txdb'` file for
#'    this purpose, or will create and save this file for future use.
#'    Similarly, the 'tx2gene' will be saved to '.tx2gene.txt'
#'    for future use.
#'    * Data are subset by 'detectedGenes' and 'detectedTx'
#'    if provided.
#'    * When both TSS and TTS data are defined, it follows rules
#'    described in the 'Details' section.
#' @param detected_genes `character` (default NULL) with optional subset
#'    of genes that are considered "detected" by some criteria.
#'    The TSS results will only include detected genes by default.
#'    See `jamses::se_detected_rows()` (Github 'jmw86069/jamses')
#'    for an example of heuristics to define detected genes in
#'    bulk RNA-seq data.
#' @param countsFromAbundance `character` string passed to
#'    `tximport::summarizeToGene()`, default 'lengthScaledTPM'.
#' @param use_signal `character` string indicating which measure of expression
#'    to use:
#'    * `"counts"` - uses data in the `counts` matrix, default.
#'    * `"abundance"` - uses data in the `abundance` matrix, TPM.
#' @param return_type `character` type of data to return:
#'    * `"tss"` (default) returns a `data.frame` with one TSS per gene
#'    for genes with non-zero observed expression.
#'    * `"bed"` as above, returns a `data.frame` with one TSS per gene
#'    for genes with non-zero observed expression, using BED columns.
#'    * `"se"` returns `SummarizedExperiment` with all TSS rows included,
#'    and `rowData()` which contains annotations for each unique TSS.
#' @param pct_thresholds `numeric` values with optional thresholds to
#'    indicate when a particular TSS met a percent-max for the given gene.
#'    * The defaults `c(50, 20)` indicate whether a TSS had expression at
#'    least 20% or 50% of the max observed expression.
#'    * The calculation is done in normal space, with exponentiated data
#'    using: `(2^x - 1)`
#' @param remove_genebody `logical` default TRUE, whether to remove
#'    transcripts recognized as full gene body (unspliced), which
#'    also have `"_gb"` suffix.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `jamba::rowGroupMeans()`
#'    for example one may customize with `useMedian=TRUE` to calculate
#'    group median expression.
#'
#' @returns Dependent upon the input arguments:
#' * When 'gtf' is provided without 'tx2tss' nor 'tx2tts',
#' or when both 'tx2tss' and 'tx2tts' are provided, both the
#' TSS and TTS sites are defined.
#' The output will be a `list` with names 'TSS' and 'TTS',
#' containing the object type as defined below.
#' * The object type is defined by argument `return_type`:
#'    * 'tss': `data.frame` default.
#'    * 'se': `SummarizedExperiment` using the input `txiTx` when also
#'    supplied as `SummarizedExperiment`. Useful with 'TxSE' produced
#'    from `platjam::import_salmon_quant()` (Github 'jmw86069/platjam')
#'    for example.
#'    * 'bed': `data.frame` in BED column format, suitable to save
#'    as BED file (without headers) or convert to `GRanges` directly, using
#'    `as(bed, "GRanges")` via the package GenomicRanges.
#'
#' @examples
#' jamsession::refresh_functions("salmontsspergene")
#' tss_df638 <- salmon_tss_per_gene(TxSE638, groups="Group", tx2tss=tx2tss, detected_genes=env638$detected_genes, verbose=TRUE)
#'
#' @export
salmon_tss_per_gene <- function(
   txiTx,
   groups = NULL,
   tx2tss = NULL,
   tx2tts = NULL,
   gtf = NULL,
   detectedGenes = NULL,
   detectedTx = NULL,
   countsFromAbundance = "lengthScaledTPM",
   use_signal = c("counts", "abundance"),
   return_type = c("tss", "bed", "se"),
   pct_thresholds = c(50, 20),
   remove_genebody = TRUE,
   verbose = FALSE,
   ...
) {
   # args
   return_type <- match.arg(return_type)
   use_signal <- match.arg(use_signal)

   # validate tss,tts,gtf input
   if (length(tx2tss) == 0 && length(tx2tts) == 0) {
      if (length(gtf) == 0) {
         stop("One must be supplied: 'tx2tss', 'tx2tt2', 'gtf'.")
      }
      tss_tts_df <- get_tss_tts_from_gtf(
         gtf = gtf,
         detectedGenes = detectedGenes,
         detectedTx = detectedTx,
         ...
      )
      tx2tss <- tss_tts_df[, c("transcript_id", "tss"), drop = FALSE]
      tx2tts <- tss_tts_df[, c("transcript_id", "tts"), drop = FALSE]
   }
   tx2tss_txid <- NULL

   #####################################
   ## Handle both TSS and TTS together
   if (length(tx2tss) > 0 && length(tx2tts) > 0) {
      if (nrow(tx2tss) != nrow(tx2tts)) {
         stop("'tx2tss' and 'tx2tts' have unequal nrow().")
      }
      if (!all(tx2tss$transcript_id == tx2tts$transcript_id)) {
         stop("'tx2tss' and 'tx2tts' have unequal rownames().")
      }
      #
      # handle tss and tts together
      tss_output <- salmon_tss_per_gene(
         txiTx = txiTx,
         groups = groups,
         tx2tss = tx2tss,
         tx2tts = NULL,
         countsFromAbundance = countsFromAbundance,
         use_signal = use_signal,
         return_type = return_type,
         pct_thresholds = pct_thresholds,
         remove_genebody = remove_genebody,
         verbose = verbose,
         ...
      )
      use_cols <- c("gene_name", "seqnames", "start", "strand")
      if (inherits(tss_output, "data.frame")) {
         use_tss <- jamba::pasteByRow(tss_output[, use_cols, drop = FALSE])
      } else if (inherits(tss_output, "SummarizedExperiment")) {
         use_tss <- jamba::pasteByRow(
            SummarizedExperiment::rowData(tss_output)[, use_cols, drop = FALSE]
         )
      } else {
         stop("'tss_output' object type not expected.")
      }
      # apply filter to TTS transcripts
      use_tts_tx <- subset(tx2tss, tss %in% use_tss)$transcript_id
      tx2tts <- subset(tx2tts, transcript_id %in% use_tts_tx)
      tts_output <- salmon_tss_per_gene(
         txiTx = txiTx,
         groups = groups,
         tx2tss = tx2tts,
         tx2tts = NULL,
         countsFromAbundance = countsFromAbundance,
         use_signal = use_signal,
         return_type = return_type,
         pct_thresholds = pct_thresholds,
         remove_genebody = remove_genebody,
         verbose = verbose,
         ...
      )
      # For now, return a list
      return(list(
         tss=tss_output,
         tts=tts_output
      ))
      #
   }

   # validate tx2tss
   tss_heading <- "tss"
   if (length(tx2tss) > 0) {
      if (!all(grepl("_chr[0-9A-Za-z_]+_[0-9]+_[-+]$", tx2tss[, 2]))) {
         stop(paste0(
            "tx2tss column 2 must have '_' delimited values: ",
            "[gene_name]_[seqnames]_[start]_[strand]"
         ))
      }
      tx2tss_txid <- tx2tss[, 1]
   } else if (length(tx2tts) > 0) {
      if (!all(grepl("_chr[0-9A-Za-z_]+_[0-9]+_[-+]$", tx2tts[, 2]))) {
         stop(paste0(
            "tx2tts column 2 must have '_' delimited values: ",
            "[gene_name]_[seqnames]_[start]_[strand]"
         ))
      }
      tss_heading <- "tts"
      if (length(tx2tss_txid) == 0) {
         tx2tss_txid <- tx2tts[, 1]
      }
   }


   # processing
   use_assay_names <- c("counts", "abundance", "length")
   TxSE638 <- NULL
   if (inherits(txiTx, "SummarizedExperiment")) {
      TxSE638 <- txiTx
      assay_names <- SummarizedExperiment::assayNames(TxSE638)
      if (!all(use_assay_names %in% assay_names)) {
         stop("assays() must contain all: counts, abundance, length")
      }
      use_tx2 <- rownames(TxSE638)
      use_tx1 <- use_tx2
      if (TRUE %in% remove_genebody) {
         use_tx1 <- jamba::unvigrep("_gb", use_tx2)
      }
      use_tx <- intersect(use_tx1, tx2tss_txid)
      if (verbose) {
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            "Recognized SummarizedExperiment input:"
         )
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            indent = 5,
            jamba::formatInt(length(use_tx2)),
            " input rows"
         )
         if (length(use_tx2) - length(use_tx1) > 0) {
            jamba::printDebug(
               "salmon_tss_per_gene(): ",
               indent = 5,
               jamba::formatInt(length(use_tx2) - length(use_tx1)),
               " genebody rows were removed"
            )
         }
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            indent = 5,
            jamba::formatInt(length(use_tx)),
            " rows found in tx2tss."
         )
         diff1 <- (length(use_tx1) - length(use_tx))
         if (diff1 > 0) {
            jamba::printDebug(
               "salmon_tss_per_gene(): ",
               indent = 5,
               jamba::formatInt(diff1),
               paste0(
                  " row",
                  ifelse(diff1 > 1, "s", ""),
                  " not found in tx2tss"
               )
            )
         }
      }

      txiTx <- lapply(
         SummarizedExperiment::assays(
            TxSE638[use_tx, ]
         )[use_assay_names],
         function(i) {
            i
         }
      )
      if (length(groups) == 0) {
         groups <- rep("Group", ncol(TxSE638))
      } else if (
         all(groups %in% colnames(SummarizedExperiment::colData(TxSE638)))
      ) {
         groups <- jamba::pasteByRow(data.frame(
            check.names = FALSE,
            SummarizedExperiment::colData(TxSE638)[, groups, drop = FALSE]
         ))
      }
   } else if (!inherits(txiTx, "list")) {
      stop("Input txiTx must be list or SummarizedExperiment.")
   }
   if (!all(use_assay_names %in% names(txiTx))) {
      stop("names(txiTx) must contain all: counts, abundance, length")
   }
   if (length(groups) == 0) {
      if (verbose) {
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            "No ",
            "groups",
            " were defined."
         )
      }
      groups <- rep("Group", ncol(txiTx[["counts"]]))
   }
   # verify signal is normal space
   if (max(txiTx$counts, na.rm = TRUE) < 50) {
      if (verbose) {
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            "Exponentiated the input counts using: ",
            "(2^x - 1)"
         )
      }
      txiTx$counts <- 2^txiTx$counts - 1
   }
   if (max(txiTx$abundance, na.rm = TRUE) < 50) {
      if (verbose) {
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            "Exponentiated the input abundance using: ",
            "(2^x - 1)"
         )
      }
      txiTx$abundance <- 2^txiTx$abundance - 1
   }

   if (length(TxSE638) == 0) {
      use_tx2 <- rownames(txiTx$counts)
      use_tx1 <- use_tx2
      if (TRUE %in% remove_genebody) {
         use_tx1 <- jamba::unvigrep("_gb", use_tx2)
      }
      use_tx <- intersect(use_tx1, tx2tss_txid)
      if (verbose) {
         jamba::printDebug("salmon_tss_per_gene(): ", "Recognized list input:")
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            indent = 5,
            jamba::formatInt(length(use_tx2)),
            " input rows"
         )
         if (length(use_tx2) - length(use_tx1) > 0) {
            jamba::printDebug(
               "salmon_tss_per_gene(): ",
               indent = 5,
               jamba::formatInt(length(use_tx2) - length(use_tx1)),
               " genebody rows were removed"
            )
         }
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            indent = 5,
            jamba::formatInt(length(use_tx)),
            " rows found in tx2tss."
         )
         diff1 <- (length(use_tx1) - length(use_tx))
         if (diff1 > 0) {
            jamba::printDebug(
               "salmon_tss_per_gene(): ",
               indent = 5,
               jamba::formatInt(diff1),
               paste0(
                  " row",
                  ifelse(diff1 > 1, "s", ""),
                  " not found in tx2tss"
               )
            )
         }
      }
   }

   if (length(use_tx) == 0) {
      stop("No rows to analyze, use verbose=TRUE for more information.")
   }
   txiTx <- lapply(txiTx[use_assay_names], function(i) {
      i[match(use_tx, rownames(i)), , drop = FALSE]
   })
   # lapply(txiTx, head, 3)

   # summarizeToTss
   txiTx$countsFromAbundance <- countsFromAbundance
   if (verbose) {
      jamba::printDebug(
         "salmon_tss_per_gene(): ",
         "Applying: ",
         "tximport::summarizeToGene()"
      )
   }
   use_colnames <- NULL
   if (length(tx2tss) > 0) {
      txiTss <- suppressMessages(
         tximport::summarizeToGene(
            txiTx,
            countsFromAbundance = countsFromAbundance,
            tx2gene = tx2tss
         )
      )
      use_colnames <- colnames(txiTss$counts)
      use_rownames <- rownames(txiTss$counts)
   }

   ## Make SummarizedExperiment for convenience
   TssSE638 <- SummarizedExperiment::SummarizedExperiment(
      assays = head(txiTss, -1),
      colData = data.frame(
         check.names = FALSE,
         row.names = use_colnames,
         Sample_ID = use_colnames,
         Group = groups
      ),
      rowData = data.frame(
         check.names = FALSE,
         gene_ref_tss = use_rownames,
         gene_name = gsub("_chr[^_]+_[^_]+_[^_]+$", "", use_rownames),
         jamba::rbindList(
            strsplit(
               gsub("^.+_(chr[^_]+_[^_]+_[^_]+$)", "\\1", use_rownames),
               "_"
            ),
            newColnames = c("seqnames", "start")
         )
      )
   )
   # log2(1 + x) transform data
   SummarizedExperiment::assays(TssSE638)$counts <- log2(
      1 +
         SummarizedExperiment::assays(TssSE638)$counts
   )
   SummarizedExperiment::assays(TssSE638)$abundance <- log2(
      1 +
         SummarizedExperiment::assays(TssSE638)$abundance
   )

   # row group mean
   if (verbose) {
      jamba::printDebug("salmon_tss_per_gene(): ",
         "Using signal: ",
         use_signal)
   }
   # useMedian=FALSE ensures that non-zero values contribute in sparse data
   rgm638 <- jamba::rowGroupMeans(
      useMedian = FALSE,
      SummarizedExperiment::assays(TssSE638)[[use_signal]],
      groups = groups,
      ...
   )
   # max group mean
   rgm638max <- apply(rgm638, 1, max)

   # summarize
   rgm638max_df <- data.frame(
      gene_name = gsub("_chr[^_]+_[^_]+_[^_]+$", "", names(rgm638max)),
      rgm638max = rgm638max
   )
   rgm638max_bygene <- jamses::shrinkDataFrame(
      rgm638max_df,
      groupBy = "gene_name",
      include_num_reps = TRUE,
      num_func = max
   )

   match_byg <- match(rgm638max_df[, 1], rgm638max_bygene[, 1])
   rgm638max_df$gene_max <- rgm638max_bygene[match_byg, 2]
   rgm638max_df$is_gene_max <- (rgm638max_df$gene_max ==
      rgm638max_df$rgm638max &
      rgm638max_df$rgm638max > 0)
   rgm638max_df$num_tsses <- rgm638max_bygene[match_byg, "num_reps"]

   ## optional summary metrics
   #
   # optional percent max expression threshold
   if (
      length(pct_thresholds) > 0 &&
         is.numeric(pct_thresholds) &&
         all(pct_thresholds >= 0 & pct_thresholds <= 1)
   ) {
      pct_thresholds <- pct_thresholds * 100
   }
   pct_thresholds <- rev(sort(unique(
      pct_thresholds[pct_thresholds > 0 & pct_thresholds < 100]
   )))
   if (length(pct_thresholds) > 0 && is.numeric(pct_thresholds)) {
      if (verbose) {
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            "Applying percent-max thresholds: ",
            pct_thresholds
         )
      }
      for (pct in pct_thresholds) {
         pct_colname <- paste0("is_", pct, "pct_gene_max")
         rgm638max_df[, pct_colname] <- (rgm638max_df$rgm638max > 0 &
            (2^rgm638max_df$rgm638max - 1) >=
               ((2^rgm638max_df$gene_max - 1) * (pct / 100)))
      }
   }

   # optional detectedGenes
   if (length(detectedGenes) > 0) {
      if (verbose) {
         jamba::printDebug(
            "salmon_tss_per_gene(): ",
            "Filtering by: ",
            "detectedGenes"
         )
      }
      rgm638max_df$is_detected_gene <- (rgm638max_df$gene_name %in%
         detectedGenes)
      rgm638max_df_tss <- subset(rgm638max_df, is_gene_max & is_detected_gene)
   } else {
      rgm638max_df_tss <- subset(rgm638max_df, is_gene_max)
   }
   if (verbose) {
      jamba::printDebug(
         "salmon_tss_per_gene(): ",
         "Identified Best TSSes: ",
         jamba::formatInt(nrow(rgm638max_df_tss))
      )
   }

   # Optionally add to SE
   if ("se" %in% return_type) {
      SummarizedExperiment::rowData(TssSE638)[, colnames(
         rgm638max_df
      )] <- (rgm638max_df)
      return(TssSE638)
   }

   use_tss_638 <- data.frame(
      check.names = FALSE,
      jamba::rbindList(
         strsplit(
            gsub("_(chr[^_]+)_([0-9]+)_([^_]+)$", "!\\1!\\2!\\3", rownames(rgm638max_df_tss)),
            "!"
         ),
         newColnames = c("gene_name", "seqnames", "start", "strand")
      )
   )
   use_tss_638$start <- as.numeric(use_tss_638$start)
   if ("bed" %in% return_type) {
      use_tss_638$end <- use_tss_638$start
      use_tss_638$start <- use_tss_638$start - 1
      use_tss_638$score <- 1
      use_tss_638_bed <- data.frame(
         check.names = FALSE,
         use_tss_638[, c("seqnames", "start", "end", "gene_name", "strand")]
      )
      return(use_tss_638_bed)
   }
   return(use_tss_638)
}

#' Get TSS and TSS sites from GTF
#' @keywords internal
#' @noRd
#' @examples
#' tss_tts_df <- get_tss_tts_from_gtf(gtf=gtf);
#'
#' tx2tss <- tss_tts_df[, c("transcript_id", "tss")]
#'
#' tx2tts <- tss_tts_df[, c("transcript_id", "tts")]
#'
get_tss_tts_from_gtf <- function(
   gtf = NULL,
   detectedGenes = NULL,
   detectedTx = NULL,
   ...
) {
   #
   if (length(gtf) == 0) {
      stop("'gtf' is required.")
   }
   refgene_txdb <- convert_gtf_to_txdb(gtf = gtf, ...)
   #
   tx2geneDF <- convert_gtf_to_tx2gene(gtf = gtf, ...)

   # transcript ranges
   txByGene <- GenomicFeatures::transcriptsBy(x = refgene_txdb, by = "gene")
   GenomicRanges::values(txByGene@unlistData)$gene_id <- rep(
      names(txByGene),
      lengths(txByGene)
   )
   genematch <- match(
      GenomicRanges::values(txByGene@unlistData)$gene_id,
      tx2geneDF[["gene_name"]]
   )
   GenomicRanges::values(txByGene@unlistData)$gene_name <- ifelse(
      is.na(genematch),
      GenomicRanges::values(txByGene@unlistData)$gene_id,
      tx2geneDF[genematch, "gene_name"]
   )
   GenomicRanges::values(txByGene@unlistData)$transcript_id <-
      GenomicRanges::values(txByGene@unlistData)$tx_name

   # Keep only flat GRanges
   txByGene <- txByGene@unlistData
   if (length(detectedGenes) > 0) {
      txByGene <- subset(
         txByGene,
         gene_name %in% detectedGenes | gene_id %in% detectedGenes
      )
   }
   if (length(detectedTx) > 0) {
      txByGene <- subset(
         txByGene,
         transcript_id %in% detectedTx
      )
   }

   # Define TSS as 1-base promoter
   tssByGene <- GenomicRanges::promoters(
      txByGene,
      upstream = 0,
      downstream = 0,
      use.names = TRUE
   )
   #
   ttsByGene <- GenomicRanges::terminators(
      txByGene,
      upstream = 0,
      downstream = 1,
      use.names = TRUE
   )

   # confirm values are still present
   GenomicRanges::values(tssByGene)$gene_name <-
      GenomicRanges::values(txByGene)$gene_name
   GenomicRanges::values(ttsByGene)$gene_name <-
      GenomicRanges::values(txByGene)$gene_name
   GenomicRanges::values(tssByGene)$transcript_id <-
      GenomicRanges::values(txByGene)$transcript_id
   GenomicRanges::values(ttsByGene)$transcript_id <-
      GenomicRanges::values(txByGene)$transcript_id

   # data.frame
   usecols <- c("gene_name", "seqnames", "start", "strand")
   tss_tts_df <- data.frame(
      transcript_id = GenomicRanges::values(txByGene)$transcript_id,
      tss = jamba::pasteByRow(
         as.data.frame(tssByGene)[, usecols, drop = FALSE]
      ),
      tts = jamba::pasteByRow(
         as.data.frame(ttsByGene)[, usecols, drop = FALSE]
      )
   )

   #
   # * `"tss"` - in the exact form `"[gene_name]_[seqnames]_[start]_[strand]"`
   # for example `"5_8S_rRNA_chr14_16057472"` is recognized to have these
   # fields: `"5_8S_rRNA", "chr14", "16057472", "+"`

   return(tss_tts_df)
}
