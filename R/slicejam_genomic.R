
#' Prepare genomic regions from GTF
#' 
#' Prepare genomic regions from GTF
#' 
#' This function takes a GTF file as input, and creates
#' a `GenomicRanges` object that represents genomic regions.
#' By default, the resulting R object is saved in a file
#' with file extension `".genome_regions.RData"`, so that
#' the data can be loaded rapidly during subsequent analyses.
#' 
#' The main components of genomic regions defined by this function:
#' 
#' * `"Promoters"`: defined with `upstream_promoter` distance upstream
#' of each gene transcript start site, and `downstream_promoter`
#' distance downstream each gene transcript start site. Note that these
#' distances are stranded, which means that negative-strand features
#' `strand="-"` will have these distances applied in the opposite
#' direction.
#' * `"TTS"`: defined as `"Transcript Termination Site"` for each
#' gene transcript. A region around the TTS is defined by applying
#' `upstream_tts` and `downstream_tts` distances, in the same manner
#' as upstream and downstream distances are applied to `"Promoters"`.
#' * `"Exons"`: defined by gene transcript exons, reduced into
#' one set of exons for each gene.
#' * `"Introns"`: defined for each gene, using gene transcript ranges
#' (from the transcript start site through the transcript termination site)
#' then subtracting exons. Introns are therefore defined as regions
#' of a gene locus which are not annotated as exons in any associated
#' gene transcripts for that gene locus.
#' 
#' Any region of the genome which is not annotated as a Promoter, TTS,
#' exon, or intron, is referred to as `"extragenic"`, although extragenic
#' regions are not defined by this function directly. Note that we favor
#' the term `"extragenic"` over `"intergenic"` because: the latter term
#' may be confused with `"intragenic"`; and may imply a region that
#' must be between two existing genes, which is not the case.
#' 
#' In more detail: The arguments `geneAttrNames`, `txAttrNames`, and
#' `geneFeatureType` are used by `splicejam::makeTx2geneFromGtf()` to
#' create a `data.frame` which relates `transcript_id` to each `gene_name`
#' by default. If the GTF file uses different attributes, adjust
#' `geneAttrNames` and `txAttrNames` to match the GTF file.
#' 
#' ### Comments on compatible files
#' 
#' The function itself is capable of reading GTF or GFF3 formatted files,
#' by calling `txdbmaker::makeTxDbFromGFF()` to do that heavy work.
#' However `splicejam::makeTx2geneFromGtf()` is required to import
#' additional gene and transcript annotations.
#' 
#' We found it is helpful for column 9 of the GTF or GFF3 file to contain
#' attributes:
#' 
#' * gene_id: a unique gene identifier
#' * gene_name: a unique gene symbol associated with each gene_id
#' * transcript_id: one or more transcript identifiers associated with
#' a gene_id
#' 
#' ### Issues with 'Name'
#' 
#' We found it problematic when the GTF or GFF3 file contains `"Name"`,
#' specifically in cases where the same value is associated with multiple
#' rows in the file. For example, a row with type "gene", `Name=APOE;`, then
#' on the next row for type "transcript", also assigns `Name=APOE;`.
#' Apparently this input causes `txdbmaker::makeTxDbFromGFF()` to assign
#' the Name value as the primary identifier for both the gene and the
#' transcript, and therefore interferes with the `tx2gene` association
#' of `gene_id` to `gene_name`, and `gene_id` to `transcript_id`.
#' 
#' This type of formatting is present in the human telomere-to-telomere
#' (T2T) v2.0 GTF file (as of November 2022).
#' 
#' In this case, we found it successful to remove the `Name` entries
#' from the GTF or GFF3 file, then call `genomic_regions_from_gtf()`,
#' potentially with `force_refresh=TRUE` to ensure the file is re-processed.
#' While inconvenient, in principle the GTF or GFF3 file is incorrect by
#' using the same `Name` value for multiple features. We did not want
#' to modify this function (yet) to accommodate an incorrect GTF or GFF3
#' input file.
#' 
#' All that said, we are incredibly grateful and happy to use
#' the human T2T v2.0 GTF gene annotations, which contain very high
#' quality gene annotations.
#' 
#' @family slicejam genome regions
#' 
#' @return `GRanges` object
#' 
#' @param gtf `character` path to a GTF file. When `gtf` is
#'    not supplied, `rdata_file` must be supplied, in which
#'    case that file is used to load pre-existing `genome_regions`
#'    GRanges object.
#' @param upstream_promoter,downstream_promoter `numeric` value
#'    which defines the distance upstream and downstream relative
#'    to each gene transcript start site, to be annotated as
#'    a `"Promoter"`.
#' @param upstream_tts,downstream_tts `numeric` value
#'    which defines the distance upstream and downstream relative
#'    to each gene transcript termination site, to be annotated as
#'    a `"TTS"`.
#' @param detectedTx `character` vector of transcripts which are
#'    "detected" and therefore used in determining the relevant
#'    gene-transcript annotations. This argument is intended to
#'    allow using a subset of transcripts, as opposed to using all
#'    annotated transcripts in the GTF file. When defined, it
#'    forces the argument `save_rdata=FALSE`.
#' @param detectedGenes `character` vector of genes which are
#'    "detected" and therefore used in determining the relevant
#'    gene-transcript annotations. This argument can be used together
#'    with `detectedTx` to define a specific subset of genes, using
#'    only the subset of detected transcripts from those genes.
#'    When defined, it forces the argument `save_rdata=FALSE`.
#' @param geneAttrNames `character` vector of attribute names to be
#'    read from the GTF file in column 9, associated at the gene level.
#'    The default values assume that `"gene_id"` and `"gene_name"`
#'    attributes are defined in the GTF file.
#' @param txAttrNames `character` vector of attribute names to be
#'    read from the GTF file in column 9, associated at the transcript
#'    level. The default values assume that `"transcript_id"`
#'    attribute is defined in the GTF file.
#' @param geneFeatureType,txFeatureType `character` vector indicating
#'    the feature types as defined in column 3 of the GTF file, associated
#'    to each gene, and transcript, respectively. These values are
#'    specifically used by `splicejam::makeTx2geneFromGtf()` to create
#'    a `data.frame` that associated transcripts to genes.
#' @param save_rdata `logical` indicating whether to save the
#'    resulting R object in a file with the file extension
#'    `".genomic_regions.RData"`.
#' @param rdata_file `character` filename, optionally used to
#'    define an input file for `genome_regions`, or a specific
#'    output file to save the RData object. When `force_refresh=TRUE`
#'    this argument is used to save new RData to a file, without
#'    re-using any pre-existing data saved in that file.
#' @param save_txdb `logical` indicating whether to save the
#'    intermediate `Txdb` R object as a SQLite database, using
#'    `AnnotationDbi::saveDb()`, using the file extension `".txdb"`.
#'    When this option is enabled, any previously stored `TxDb`
#'    will be re-used, unless `force_refresh=TRUE`.
#' @param save_bed `logical` indicating whether to save the
#'    genome_regions also in BED format.
#' @param force_refresh `logical` indicating whether to force
#'    a full refresh of the processing steps in this function.
#'    When `force_refresh=TRUE`, the `rdata_file` input file
#'    is not re-used, but is created and will overwrite
#'    the `rdata_file` if it exists.
#'    When `force_refresh=TRUE`: the GTF file will be read;
#'    a new TxDb object will be created, optionally saved if
#'    `save_txdb=TRUE`; and new RData file will be saved if
#'    `save_rdata=TRUE`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @export
genomic_regions_from_gtf <- function
(gtf=NULL,
 upstream_promoter=2000,
 downstream_promoter=200,
 upstream_tts=1000,
 downstream_tts=1000,
 detectedTx=NULL,
 detectedGenes=NULL,
 geneAttrNames=c("gene_id", "gene_name"),
 txAttrNames=c("transcript_id"),
 geneFeatureType="exon",
 txFeatureType="exon",
 save_rdata=TRUE,
 rdata_file=NULL,
 save_txdb=TRUE,
 save_bed=TRUE,
 force_refresh=FALSE,
 verbose=FALSE,
 ...)
{
   ##
   ## genome_regions may be saved to a file already
   if (length(gtf) == 0 || !file.exists(gtf)) {
      if (length(rdata_file) == 0 || !file.exists(rdata_file)) {
         stop("Must supply either a GTF file, or a RData file with 'genome_regions' object.");
      }
      if (force_refresh) {
         stop("Cannot force_refresh when GTF is not supplied.");
      }
      if (length(detectedTx) > 0 || length(detectedGenes) > 0) {
         warning("Note detectedTx and detectedGenes are ignored when no GTF is supplied.");
      }
   }
   if (length(gtf) > 0) {
      if (!file.exists(gtf)) {
         stop(paste0("GTF file not found:", gtf));
      }
      gtf_base <- gsub("[.](gff|gtf|gff3)([-_.](zip|gz|Z|tgz|tar.gz)|)$",
         "",
         ignore.case=TRUE,
         gtf);
      if (length(rdata_file) == 0) {
         ## Create suitable filename for rdata_file
         short_dist_label <- function(x){
            if (x >= 1000) {
               paste0(x/1000, "kb")
            } else {
               paste0(x, "b")
            }
         }
         if (exists("mask_region") && length(mask_regions) > 0) {
            mask_ext <- "_mask";
         } else {
            mask_ext <- "";
         }
         save_ext <- paste0(
            ".",
            short_dist_label(upstream_promoter),
            "_tss_",
            short_dist_label(downstream_promoter),
            ".",
            short_dist_label(upstream_tts),
            "_tts_",
            short_dist_label(downstream_tts),
            mask_ext,
            ".genome_regions.RData");
         gtf_gr_file <- paste0(gtf_base, save_ext);
         if (!file.exists(gtf_gr_file) && file.exists(basename(gtf_gr_file))) {
            gtf_gr_file <- basename(gtf_gr_file);
         }
         if (length(detectedTx) > 0 || length(detectedGenes) > 0) {
            save_rdata <- FALSE;
         }
      } else {
         gtf_gr_file <- rdata_file;
      }
   }
   if (file.exists(gtf_gr_file)) {
      if (force_refresh) {
         if (verbose) {
            jamba::printDebug("genomic_regions_from_gtf(): ",
               "Refreshing genome_regions into existing RData file:",
               gtf_gr_file);
         }
      } else if (file.exists(gtf_gr_file) && save_rdata) {
         if (verbose) {
            jamba::printDebug("genomic_regions_from_gtf(): ",
               "Reloading genome_regions from RData file:",
               gtf_gr_file);
         }
         genome_regions_o <- load(gtf_gr_file);
         #if (!exists("genome_regions")) {
         if (!"genome_regions" %in% genome_regions_o) {
            jamba::printDebug("The RData file '",
               gtf_gr_file,
               "' did not contain R object ",
               '"genome_regions"', ".");
            jamba::printDebug("Instead, it contained these R objects:",
               genome_regions_o);
            stop("Invalid RData file, no genome_regions object was found.");
         }
         if (!"rdata_file" %in% names(attributes(genome_regions))) {
            attr(genome_regions, "rdata_file") <- gtf_gr_file;
         }
         # only re-use genome_regions if we do not need to subset the data
         if (length(detectedTx) == 0 && length(detectedGenes) == 0) {
            return(genome_regions);
         }
      }
   }

   ## If genome_regions is not define, create it
   # if we need to create genome_regions, we need the txdb
   #if (!exists("genome_regions")) {
   txdb_file <- paste0(gtf_base, ".txdb");
   refgene_txdb <- NULL;
   if (!force_refresh && save_txdb) {
      if (!file.exists(txdb_file) && file.exists(basename(txdb_file))) {
         txdb_file <- basename(txdb_file);
      }
      if (file.exists(txdb_file)) {
         if (verbose) {
            jamba::printDebug("genomic_regions_from_gtf(): ",
               "Loading existing txdb file:",
               txdb_file);
         }
         refgene_txdb <- AnnotationDbi::loadDb(txdb_file);
      }
   }
   if (length(refgene_txdb) == 0) {
      ## 15-20 seconds for human GTF
      if (verbose) {
         jamba::printDebug("genomic_regions_from_gtf(): ",
            "Creating txdb from gtf:",
            gtf);
      }
      refgene_txdb <- jamba::call_fn_ellipsis(
         txdbmaker::makeTxDbFromGFF,
         file=gtf,
         ...);
      if (save_txdb) {
         tryCatch({
            AnnotationDbi::saveDb(refgene_txdb,
               file=txdb_file);
         }, error=function(e){
            txdb_file <- basename(txdb_file);
            AnnotationDbi::saveDb(refgene_txdb,
               file=txdb_file);
         });
         if (verbose) {
            if (file.exists(txdb_file)) {
               jamba::printDebug("genomic_regions_from_gtf(): ",
                  "Saved txdb to file:",
                  txdb_file);
            } else {
               jamba::printDebug("genomic_regions_from_gtf(): ",
                  "Unable to save txdb to file:",
                  txdb_file);
            }
         }
      }
   }
   #}
   
   ## tx2geneDF
   if (!exists("tx2geneDF")) {
      if (verbose) {
         jamba::printDebug("genomic_regions_from_gtf(): ",
            "Creating tx2geneDF");
      }
      tx2gene_file <- paste0(
         gsub("[.](gff|gff3|gtf)(|[.]gz)$", "", ignore.case=TRUE, gtf),
         ".tx2gene.txt");
      if (!file.exists(tx2gene_file) && file.exists(basename(tx2gene_file))) {
         tx2gene_file <- basename(tx2gene_file);
      }
      if (file.exists(tx2gene_file)) {
         tx2geneDF <- data.table::fread(tx2gene_file,
            sep="\t",
            data.table=FALSE);
      } else {
         tx2geneDF <- splicejam::makeTx2geneFromGtf(gtf,
            geneAttrNames=geneAttrNames,
            txAttrNames=txAttrNames,
            geneFeatureType=geneFeatureType,
            verbose=verbose,
            txFeatureType=txFeatureType);
         # save to a file
         tryCatch({
            data.table::fwrite(tx2geneDF,
               file=tx2gene_file,
               sep="\t");
         }, error=function(e){
            data.table::fwrite(tx2geneDF,
               file=basename(tx2gene_file),
               sep="\t");
         });
      }
   }
   gene_colname <- head(intersect(geneAttrNames, colnames(tx2geneDF)), 1);
   if (length(gene_colname) == 0) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "geneAttrNames were not present in tx2geneDF.");
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "input geneAttrNames:",
         geneAttrNames);
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "colnames(tx2geneDF):",
         colnames(tx2geneDF));
   }
   geneAttrNames <- intersect(geneAttrNames, colnames(tx2geneDF));
   tx_colname <- head(intersect(txAttrNames, colnames(tx2geneDF)), 1);
   if (length(tx_colname) == 0) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "txAttrNames were not present in tx2geneDF.");
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "input txAttrNames:",
         txAttrNames);
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "colnames(tx2geneDF):",
         colnames(tx2geneDF));
   }
   txAttrNames <- intersect(txAttrNames, colnames(tx2geneDF));
   
   ## Process detectedTx and detectedGenes
   if (length(detectedTx) > 0) {
      if (verbose) {
         jamba::printDebug("genomic_regions_from_gtf(): ",
            "Processing the supplied detectedTx");
      }
      detectedTx <- intersect(detectedTx,
         tx2geneDF[[tx_colname]]);
      tx2geneDF <- subset(tx2geneDF, tx2geneDF[[tx_colname]] %in% detectedTx);
   }
   if (length(detectedGenes) > 0) {
      if (verbose) {
         jamba::printDebug("genomic_regions_from_gtf(): ",
            "Processing the supplied detectedGenes");
      }
      detectedGenesL <- lapply(geneAttrNames, function(gene_attr){
         (tx2geneDF[[gene_attr]] %in% detectedGenes)
      })
      detectedGenesT <- Reduce("|", detectedGenesL);
      if (!any(detectedGenesT)) {
         jamba::printDebug("detectedGenes were not recognized in geneAttrNames:",
            geneAttrNames);
         jamba::printDebug("head(detectedGenes):",
            head(detectedGenes));
         jamba::printDebug("head(tx2geneDF):");
         print(head(tx2geneDF));
         stop("detectedGenes were not found in tx2geneDF.");
      }
      tx2geneDF <- subset(tx2geneDF, detectedGenesT);
      detectedTx <- unique(tx2geneDF[[tx_colname]]);
   }
   if (verbose) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "nrow(tx2geneDF):", jamba::formatInt(nrow(tx2geneDF)));
   }
   
   ## Exons
   if (length(detectedTx) > 0) {
      ## assemble exons only for detected transcripts
      if (verbose) {
         jamba::printDebug("genomic_regions_from_gtf(): ",
            "Assembling exons for detectedTx gene transcripts.");
      }
      ## 3-5 seconds
      exonsByTx <- GenomicFeatures::exonsBy(refgene_txdb,
         by="tx",
         use.names=TRUE);
      # subset for entries contained in tx2geneDF
      ikeep <- (names(exonsByTx) %in% tx2geneDF[[tx_colname]]);
      if (!all(ikeep)) {
         exonsByTx <- exonsByTx[ikeep];
      }
      tx_match <- match(names(exonsByTx),
         tx2geneDF[[tx_colname]]);
      GenomicRanges::values(exonsByTx@unlistData)[[gene_colname]] <- rep(
         tx2geneDF[tx_match, gene_colname],
         IRanges::elementNROWS(exonsByTx));
      exonsByGene <- GenomicRanges::reduce(
         GenomicRanges::split(exonsByTx@unlistData,
            GenomicRanges::values(exonsByTx@unlistData)[[gene_colname]]));
      gene_match <- match(names(exonsByGene),
         tx2geneDF[[gene_colname]]);
      for (gene_attr in geneAttrNames) {
         GenomicRanges::values(exonsByGene@unlistData)[[gene_attr]] <- rep(
            tx2geneDF[gene_match, gene_attr],
            IRanges::elementNROWS(exonsByGene))
      }
   } else {
      ## assemble exons for all transcripts
      if (verbose) {
         jamba::printDebug("genomic_regions_from_gtf(): ",
            "Assembling exons for all gene transcripts.");
      }
      exonsByGene <- GenomicFeatures::exonsBy(refgene_txdb,
         by="gene");
      # 28nov2022: find matching gene colname
      for (geneAttrName in intersect(geneAttrNames, colnames(tx2geneDF))) {
         if (all(names(exonsByGene) %in% tx2geneDF[[geneAttrName]])) {
            gene_colname <- geneAttrName;
         }
      }
      exon_match <- match(names(exonsByGene),
         tx2geneDF[[gene_colname]]);
      ## define gene_name and gene_id to each exon GRanges
      for (gene_attr in geneAttrNames) {
         GenomicRanges::values(exonsByGene@unlistData)[[gene_attr]] <- rep(
            tx2geneDF[exon_match, gene_attr],
            IRanges::elementNROWS(exonsByGene));
      }
   }
   if (verbose) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "head(exonsByGene):");
      print(head(exonsByGene));
   }
   
   ## transcript ranges
   if (verbose) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "Creating txByGene.");
   }
   txByGene <- GenomicFeatures::transcriptsBy(refgene_txdb);
   
   # 28nov2022
   # ensure tx_name values are unique...
   GenomicRanges::values(txByGene@unlistData)[,"tx_name"] <- jamba::makeNames(
      GenomicRanges::values(txByGene@unlistData)[,"tx_name"])
   
   # rename colnames
   GenomicRanges::values(txByGene@unlistData) <- jamba::renameColumn(
      GenomicRanges::values(txByGene@unlistData),
      from=c("tx_id", "tx_name"),
      to=c("internal_tx_id", tx_colname));
   if (verbose) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "head(txByGene):");
      print(head(txByGene));
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "head(tx2geneDF):");
      print(head(tx2geneDF));
   }
   # subset for tx
   ikeep <- (GenomicRanges::values(txByGene@unlistData)[[tx_colname]] %in% tx2geneDF[[tx_colname]]);
   if (!all(ikeep)) {
      if (verbose) {
         jamba::printDebug("genomic_regions_from_gtf(): ",
            "Subsetting txByGene for detectedTx, table(ikeep):");
      }
      txByTx <- subset(txByGene@unlistData, ikeep);
      tx_match <- match(GenomicRanges::values(txByTx)[[tx_colname]],
         tx2geneDF[[tx_colname]]);
      for (gene_attr in geneAttrNames) {
         GenomicRanges::values(txByTx)[[gene_attr]] <- tx2geneDF[tx_match, gene_attr];
      }
      txByGene <- GenomicRanges::split(txByTx,
         GenomicRanges::values(txByTx)[[gene_colname]]);
   } else {
      tx_match <- match(GenomicRanges::values(txByGene@unlistData)[[tx_colname]],
         tx2geneDF[[tx_colname]]);
      for (gene_attr in geneAttrNames) {
         GenomicRanges::values(txByGene@unlistData)[[gene_attr]] <- tx2geneDF[tx_match, gene_attr];
      }
   }

   ## TTS per transcript range, extend -1000,+1000 around TTS
   if (verbose) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "Creating ttsByTx.");
   }
   ttsByTx <- GenomicRanges::flank(
      txByGene@unlistData,
      start=FALSE,
      width=-1);
   ttsByTx <- GenomicRanges::punion(
      GenomicRanges::flank(ttsByTx,
         width=-upstream_tts,
         both=FALSE,
         start=FALSE),
      GenomicRanges::flank(ttsByTx,
         width=downstream_tts,
         both=FALSE,
         start=FALSE)
   )
   GenomicRanges::values(ttsByTx) <- GenomicRanges::values(txByGene@unlistData);
   # tts reduce() per gene
   ttsByGene <- GenomicRanges::split(ttsByTx,
      GenomicRanges::values(ttsByTx)[[gene_colname]]);
   ttsByGeneRed <- GenomicRanges::reduce(ttsByGene);
   GenomicRanges::values(ttsByGeneRed@unlistData)[[gene_colname]] <- rep(
      names(ttsByGeneRed),
      lengths(ttsByGeneRed));
   gene_match <- match(GenomicRanges::values(ttsByGeneRed@unlistData)[[gene_colname]],
      GenomicRanges::values(ttsByGene@unlistData)[[gene_colname]])
   for (gene_attr in geneAttrNames) {
      GenomicRanges::values(ttsByGeneRed@unlistData)[[gene_attr]] <- GenomicRanges::values(
         ttsByGene@unlistData)[[gene_attr]][gene_match];
   }
   ttsByTxRed <- ttsByGeneRed@unlistData;
   
   ## promoters
   ## Use default values for GenomicFeatures::promoters()
   ## but store them here to make sure they do not change
   if (verbose) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "Creating promoters_gr.");
   }
   # use flank() instead of promoters() to keep same annotation
   # in the same order as txByGene@unlistData
   promoters_gr <- GenomicRanges::flank(
      txByGene@unlistData,
      start=TRUE,
      width=-1);
   promoters_gr <- GenomicRanges::punion(
      GenomicRanges::flank(promoters_gr,
         width=-upstream_promoter,
         both=FALSE,
         start=TRUE),
      GenomicRanges::flank(promoters_gr,
         width=downstream_promoter,
         both=FALSE,
         start=TRUE)
   )
   GenomicRanges::values(promoters_gr) <- GenomicRanges::values(txByGene@unlistData);
   # promoters_gr shrink() per gene
   # promoters_gr shrink() per gene
   gene_chr_start_end <- paste0(
      GenomicRanges::values(promoters_gr)[[gene_colname]], "_",
      GenomicRanges::seqnames(promoters_gr), "_",
      GenomicRanges::start(promoters_gr), "_",
      GenomicRanges::end(promoters_gr))
   if (any(duplicated(gene_chr_start_end))) {
      promoters_gr_shr <- promoters_gr[match(unique(gene_chr_start_end), gene_chr_start_end)];
   } else {
      promoters_gr_shr <- promoters_gr
   }
   # promoters_gr reduce() per gene
   if (1 == 2) {
      # this section is disabled for now but kept for future use
      tssByGene <- GenomicRanges::split(promoters_gr,
         GenomicRanges::values(promoters_gr)[[gene_colname]]);
      tssByGeneRed <- GenomicRanges::reduce(tssByGene);
      GenomicRanges::values(tssByGeneRed@unlistData)[[gene_colname]] <- rep(
         names(tssByGeneRed),
         lengths(tssByGeneRed));
      gene_match <- match(GenomicRanges::values(tssByGeneRed@unlistData)[[gene_colname]],
         GenomicRanges::values(tssByGene@unlistData)[[gene_colname]])
      for (gene_attr in geneAttrNames) {
         GenomicRanges::values(tssByGeneRed@unlistData)[[gene_attr]] <- GenomicRanges::values(
            tssByGene@unlistData)[[gene_attr]][gene_match];
      }
      promoters_gr_red <- tssByGeneRed@unlistData;
   }
   
   ## Introns
   ## defined everything within a transcript range that is not an exon
   if (verbose) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "Creating intronsByGene.");
   }
   tx_match <- match(names(txByGene), names(exonsByGene));
   intronsByGene <- GenomicRanges::setdiff(txByGene,
      exonsByGene[tx_match]);
   intron_match <- match(names(intronsByGene),
      tx2geneDF[[gene_colname]]);
   for (gene_attr in geneAttrNames) {
      GenomicRanges::values(intronsByGene@unlistData)[[gene_attr]] <- rep(
         tx2geneDF[intron_match, gene_attr],
         IRanges::elementNROWS(intronsByGene));
   }

   ## Assemble each layer of genomic region
   promoter_name <- paste0("Promoters (-",
      upstream_promoter,
      ",+",
      downstream_promoter,
      ")");
   tts_name <- paste0("TTS (-",
      upstream_tts,
      ",+",
      downstream_tts,
      ")");
   
   ## define genomic_regions
   if (verbose) {
      jamba::printDebug("genomic_regions_from_gtf(): ",
         "Creating genome_regions");
   }
   genome_regions_l <- list(
      promoters=promoters_gr_shr[,c("gene_name", "gene_id")],
      exons=exonsByGene@unlistData[,c("gene_name", "gene_id")],
      introns=intronsByGene@unlistData[,c("gene_name", "gene_id")],
      tts=ttsByTxRed[,c("gene_name", "gene_id")]);
   names(genome_regions_l)[1] <- promoter_name;
   names(genome_regions_l)[4] <- tts_name;

   ## combine list elements into one GRanges object
   genome_regions <- GenomicRanges::GRangesList(genome_regions_l)@unlistData;
   GenomicRanges::values(genome_regions)$feature_type <- rep(names(genome_regions_l),
      lengths(genome_regions_l));

   ## Save to the RData file
   if (save_rdata && length(detectedTx) == 0) {
      gtf_gr_file <- tryCatch({
         save(list=c("genome_regions", "tx2geneDF"),
            file=gtf_gr_file);
         if (verbose) {
            jamba::printDebug("genomic_regions_from_gtf(): ",
               "Genome regions were saved for later re-use:",
               gtf_gr_file);
         }
         gtf_gr_file
      }, error=function(e){
         save(list=c("genome_regions", "tx2geneDF"),
            file=basename(gtf_gr_file));
         gtf_gr_file <- basename(gtf_gr_file);
         if (verbose) {
            jamba::printDebug("genomic_regions_from_gtf(): ",
               "Genome regions were saved for later re-use:",
               gtf_gr_file);
         }
         gtf_gr_file
      });
      attr(genome_regions, "rdata_file") <- gtf_gr_file;
   }
   bed_file <- NULL;
   if (length(save_bed) > 0) {
      if (is.logical(save_bed) && save_bed && length(detectedTx) == 0) {
         bed_file <- paste0(
            gsub("[.]rdata$",
               "",
               ignore.case=TRUE,
               gtf_gr_file),
            ".bed");
      } else if (is.character(save_bed)) {
         bed_file <- save_bed;
      }
   }
   if (length(bed_file) > 0) {
      if (!file.exists(bed_file) || force_refresh) {
         attr_colnames <- c(
            head(
               jamba::provigrep(c("gene.*name", "."),
                  c(geneAttrNames, txAttrNames)), 
               1),
            "feature_type");
         gr_id <- paste0(GenomicRanges::seqnames(genome_regions), ":",
            GenomicRanges::start(genome_regions), "-",
            GenomicRanges::end(genome_regions));
         genome_regions_names <- jamba::pasteByRow(GenomicRanges::values(
            genome_regions)[,attr_colnames],
            sep="|");
         names(genome_regions) <- genome_regions_names;
         rtracklayer::export.bed(
            object=GenomicRanges::sort(genome_regions, ignore.strand=TRUE),
            con=bed_file);
         if (verbose) {
            jamba::printDebug("genomic_regions_from_gtf(): ",
               "Genome regions were saved in BED format:",
               bed_file);
         }
      }
   }
   return(genome_regions);
}

#' Create BED format from stats results data.frame
#' 
#' Create BED format from stats results data.frame
#' 
#' This function is specific to the slicejam package, intended
#' to convert a stats result `data.frame` into a `BED` format.
#' 
#' @param statsdf `data.frame` that contains results
#'    following a statistical contrast.
#' @param MGM `numeric` value indicating the max group mean
#'    threshold applied to the analyzed data.
#' @param `character` value indicating results to save: `"all"`
#'    will save results for all rows in `statsdf`; `"hits"`
#'    will save only those results that are considered
#'    statistical hits using the thresholds applied.
#' @param ... additional arguments are ignored.
#' 
#' @export
statsdf2bed <- function
(statsdf,
 MGM=NULL,
 type=c("all", "hits"),
 ...)
{
   ## This function converts Stats data.frame to BED format
   ## using the BED name field to encode additional info
   type <- match.arg(type);
   
   hit_cols <- jamba::vigrep("^hit ", colnames(statsdf));
   hit_cols_names <- gsub("^.+ ", "", hit_cols);
   hit_cols_names2 <- gsub("_v2", "",
      jamba::makeNames(hit_cols_names));
   names(hit_cols) <- hit_cols_names2;
   hit_cols_use <- hit_cols[jamba::unvigrep("_v[0-9]+$", hit_cols_names2)];
   
   ## Optionally subset for stats hits
   if ("hits" %in% type) {
      hits_l <- lapply(hit_cols_use, function(i){statsdf[[i]] != 0})
      hits_v <- Reduce("|", hits_l);
      statsdf <- subset(statsdf, hits_v);
   }
   if (nrow(statsdf) == 0) {
      return(NULL)
   }
   
   hits_im <- as.matrix(statsdf[,hit_cols_use,drop=FALSE]);
   hits_score <- rowSums(abs(hits_im));
   ## Make BED name
   hits_vl <- lapply(jamba::nameVectorN(hit_cols_use), function(i){
      j <- statsdf[[hit_cols[i]]];
      ifelse(j == 0,
         "",
         paste0(j, "_", i))
   });
   ## _annoTSS_gene_name
   mgmpg <- paste0("mgm", MGM, "_annoTSS_gene_name");
   mgmpgid <- paste0("mgm", MGM, "_annoTSS_gene_id");
   mgmpd <- paste0("mgm", MGM, "_annoTSS_distance");
   
   statsdf$ft_winner <- gsub("[ ()].*$",
      "", 
      statsdf$feature_type_winner);
   
   hits_vm <- do.call(cbind, c(
      lapply(jamba::nameVector(c("Gene", "ft_winner", mgmpg, mgmpd)), function(i){
         rmNA(naValue="", statsdf[[i]])
      }),
      hits_vl
   ));
   hits_label <- jamba::pasteByRow(hits_vm, sep="|");
   statsdf$label <- hits_label;
   statsdf$score <- hits_score;
   statsdf$strand <- "+";
   statsbed <- statsdf[,c("Chr", "Start", "End", "label", "score", "strand"),drop=FALSE];
   return(statsbed);
}

#' Annotate GRanges by genome_regions
#' 
#' Annotate GRanges by genome_regions
#' 
#' This function uses the `genome_regions` data defined by
#' `genomic_regions_from_gtf()` to annotate `GRanges` object `gr`.
#' 
#' It performs three levels of annotation:
#' 
#' 1. Direct overlap. Any overlapping region in `genome_regions`
#' is added as an annotation column, where multiple regions
#' are concatenated by commas.
#' 2. "Winner" overlap. When there are multiple overlapping
#' regions from step 1, the annotation(s) from the best `feature_type`
#' are called "winner" and appended with column suffix `"_winner"`.
#' 3. Nearest gene. The annotation of any overlapping gene body,
#' or nearest gene body to each feature in `gr`. Columns
#' have the prefix `"nearest_"`, and the distance is stored as
#' `"nearest_gene_distance"`.
#' 
#' When `mask_regions` is supplied, one additional column is
#' added `"mask_region"` with either `TRUE` or `FALSE`.
#' 
#' @import data.table
#' 
#' @family slicejam genome regions
#' 
#' @param gr `GRanges` object to be annotated. The `values(gr)`
#'    should contain a colname that matches `name_colname`,
#'    otherwise names will be created for each entry in
#'    colname `"name"`, see `name_colname` below.
#' @param genome_regions `GRanges` object as produced by
#'    `genomic_regions_from_gtf()`. The `values(genome_regions)`
#'    should contain a colname that matches `gene_name_colname`,
#'    and the feature type should be stored in a colname
#'    `feature_type_colname`.
#' @param name_colname `character` value that matches one colname
#'    in `values(gr)`. When no colname is supplied, a colname `"name"`
#'    is created with dummy values with the format `"gr00001"`.
#' @param gene_name_colname `character` value that matches one colname
#'    in `values(genome_regions)`.
#' @param feature_type_colname `character` value that matches one colname
#'    in `values(genome_regions)`, and which contains feature types.
#' @param gene_id_colname `character` optional value of one colname in
#'    `values(genome_regions)` to be included alongside gene annotations.
#' @param mask_regions one of the following input formats:
#'    * `character` vector with file or files that contain mask regions
#'    in BED format, only the regions are retained without further annotation
#'    * `GRanges` object containing mask regions. Currently the name is ignored.
#'    * `GRangesList` object, which is used as `mask_regions@unlistData` to
#'    convert that format to `GRanges` for internal use. No other annotations
#'    are used, all regions are considered `"mask"`.
#' @param feature_grep_order `character` vector of grep patterns used
#'    to match values in `feature_type_colname`, to define the priority
#'    of feature types to use for the "winner".
#' @param include_type `character` vector indicating which of the three
#'    annotation phases to include:
#'    1. `"overlap"` annotates each region by direct overlap with `genome_regions`
#'    2. `"winner"` annotates each region by direct overlap, using only
#'    the feature_type_winner for each entry in `gr`.
#'    3. `"nearest_gene"` annotates each entry in `gr` by the nearest gene
#'    in `genome_regions` alongside the distance to nearest gene.
#' @param verbose `logical` indicating whether to print verbose output.
#'    More detailed output is printed when `verbose=2`.
#' @param ... additional arguments are ignored.
#' 
#' @export
annotate_gr_by_genome_region <- function
(gr,
 genome_regions,
 name_colname="name",
 gene_name_colname="gene_name",
 feature_type_colname="feature_type",
 gene_id_colname="gene_id",
 mask_regions=NULL,
 feature_grep_order=c("promoter", 
    "tts",
    "exon",
    "intron",
    "extragenic",
    "intergenic",
    "."),
 include_type=c("overlap",
    "winner",
    "nearest_gene"),
 sort_optimization=c("fast",
    "global"),
 verbose=FALSE,
 ...)
{
   # validate input
   include_type <- match.arg(include_type,
      several.ok=TRUE);
   sort_optimization <- match.arg(sort_optimization);
   
   name_colname <- intersect(name_colname,
      colnames(GenomicRanges::values(gr)));
   if (length(name_colname) == 0) {
      name_colname <- "name";
      GenomicRanges::values(gr)[[name_colname]] <- paste0("gr",
         jamba::padInteger(seq_along(gr)));
   }
   
   gene_name_colname <- head(intersect(gene_name_colname,
      colnames(GenomicRanges::values(genome_regions))), 1);
   feature_type_colname <- head(intersect(feature_type_colname,
      colnames(GenomicRanges::values(genome_regions))), 1);
   if (length(gene_name_colname) == 0 || length(feature_type_colname) == 0) {
      stop("gene_name_colname and feature_type_colname must be present in colnames(values(genome_regions))");
   }
   gene_id_colname <- intersect(gene_id_colname,
      colnames(GenomicRanges::values(genome_regions)));
   
   ## Expand genome_regions if there are multi-gene features
   # begin to annotate peaks by genome_regions
   
   ## Expand comma-delimited gene_name values if present in genome_regions
   if (jamba::igrepHas(",", GenomicRanges::values(genome_regions)[[gene_name_colname]])) {
      gr_expl <- strsplit(
         GenomicRanges::values(genome_regions)[[gene_name_colname]], ",");
      gr_expi <- rep(seq_along(genome_regions),
         lengths(gr_expl));
      genome_regions_exp <- unname(genome_regions)[gr_expi];
      GenomicRanges::values(genome_regions_exp)[[gene_name_colname]] <- unname(unlist(gr_expl));
      
      ## Remove spaces and commas from feature_type column
      genome_regions_exp <- genome_regions;
      GenomicRanges::values(genome_regions_exp)[[feature_type_colname]] <- gsub("[, ]+",
         "",
         GenomicRanges::values(genome_regions_exp)[[feature_type_colname]]);
   } else {
      genome_regions_exp <- genome_regions;
   }

   ## combine gene_name with feature_type
   GenomicRanges::values(genome_regions_exp)[,"gene_feature_type"] <- pasteByRow(
      GenomicRanges::values(genome_regions_exp)[,c(gene_name_colname, feature_type_colname)],
      sep=" ");

   ##################################################
   ## Annotate peaks by overlapping region
   if (any(c("overlap", "winner") %in% include_type)) {
      if (verbose) {
         jamba::printDebug("annotate_gr_by_genome_region(): ",
            "Annotate peaks by overlapping region");
      }
      fco <- GenomicRanges::findOverlaps(gr,
         ignore.strand=TRUE,
         genome_regions_exp);
      qh1 <- S4Vectors::queryHits(fco);
      sh1 <- S4Vectors::subjectHits(fco);
      grdt <- data.table::data.table(
         fc=GenomicRanges::values(gr[qh1])[[name_colname]],
         gene_name=GenomicRanges::values(
            genome_regions_exp[sh1])[[gene_name_colname]],
         feature_type=factor(
            GenomicRanges::values(
               genome_regions_exp[sh1])[[feature_type_colname]],
            levels=jamba::provigrep(feature_grep_order,
               unique(GenomicRanges::values(
                  genome_regions_exp[sh1])[[feature_type_colname]]))),
         gene_feature_type=GenomicRanges::values(
            genome_regions_exp[sh1])$gene_feature_type);
      # optionally rename columns - not necessary for now
      if (FALSE) {
         grdt <- jamba::renameColumn(grdt,
            from=c("gene_name",
               "feature_type"),
            to=c(gene_name_colname,
               feature_type_colname));
      }

      # add gene_id if defined
      if (length(gene_id_colname) == 1) {
         grdt$gene_id <- GenomicRanges::values(
            genome_regions_exp[sh1])[[gene_id_colname]];
      }
      grd_colnames <- S4Vectors::intersect(
         c("gene_name",
            "gene_id",
            "feature_type",
            "gene_feature_type"),
         colnames(grdt));
      if (verbose) {
         jamba::printDebug("annotate_gr_by_genome_region(): ",
            "head(grdt):");
         print(head(grdt));
      }
      
      ## Iterate each column, combine multi-row features into one row
      ## then annotate using unique, sorted, comma-delimited values
      if ("overlap" %in% include_type) {
         if (verbose) {
            jamba::printDebug("annotate_gr_by_genome_region(): ",
               "grd_colnames:", grd_colnames);
         }
         grd_vals <- lapply(jamba::nameVector(grd_colnames), function(i){
            icols <- unname(c("fc", i));
            if (verbose) {
               jamba::printDebug("annotate_gr_by_genome_region(): ",
                  "icols:", icols);
               jamba::printDebug("annotate_gr_by_genome_region(): ",
                  "colnames(grdt):", colnames(grdt));
            }
            ## alternate syntax without "..icols"
            # idt <- grdt[, ..icols];
            idt <- grdt[, icols, with=FALSE];
            ## expand comma-delimited entries
            grdt0ft <- unique(jamba::mixedSortDF(idt, byCols=1:2)[, c(1, 2)]);
            grv <- S4Vectors::unstrsplit(
               IRanges::CharacterList(
                  split(as.character(grdt0ft[[2]]), grdt0ft[[1]])),
               sep=",");
            grv;
         });
         ## Add each column to gr
         for (i in names(grd_vals)) {
            imatch <- match(
               names(grd_vals[[i]]),
               GenomicRanges::values(gr)[[name_colname]]);
            if ("feature_type" == i) {
               GenomicRanges::values(gr)[, i] <- unname(
                  c("extragenic", grd_vals[[i]][1])[1]);
            } else {
               GenomicRanges::values(gr)[, i] <- unname(
                  c(NA, grd_vals[[i]][1])[1]);
            }
            GenomicRanges::values(gr)[imatch, i] <- grd_vals[[i]];
         }
      }
   }
   
   ## is_duplicate()
   is_duplicate <- function(x) {
      duplicated(x) | duplicated(x,
         fromLast=TRUE);
   }
   
   ##################################################
   ## sort for best feature_type
   if ("winner" %in% include_type) {
      if (verbose) {
         jamba::printDebug("annotate_gr_by_genome_region(): ",
            "Annotate winner");
      }
      
      # sort by feature_type using one of two methods
      grdt_fac <- factor(grdt$feature_type,
         levels=jamba::provigrep(feature_grep_order,
            unique(grdt$feature_type)));
      if (FALSE) {# && "global" %in% sort_optimization) {
         grdt$feature_type <- grdt_fac;
         sort_time <- system.time({
            grdt0 <- jamba::mixedSortDF(grdt,
               byCols=c("fc",
                  "feature_type",
                  "gene_name",
                  "gene_id"
                  ));
         }, gcFirst=FALSE);
      } else {#if ("fast" %in% sort_optimization) {
         sort_time <- system.time({
            grdt0 <- grdt[order(grdt_fac)];
         }, gcFirst=FALSE);
      }
      if (verbose) {
         jamba::printDebug("",
            "elapsed sort_time (sec): ", sort_time[["elapsed"]]);
      }
      
      ## find first row per peak
      if (TRUE) {
         grdt0_bestft <- jamba::nameVector(subset(grdt0,
            !duplicated(fc))[, c("feature_type", "fc"), drop=FALSE]);
      } else {
         grdt0_peak <- grdt0[["fc"]];
         grdt0_peaku <- unique(grdt0_peak);
         whichu <- match(grdt0_peaku, grdt0_peak);
         ## First row per peak
         grdt0_pk <- grdt0[whichu,][["fc"]];
         grdt0_bestft <- grdt0[whichu,][["feature_type"]];
         names(grdt0_bestft) <- grdt0_pk;
      }
      
      ## subset for peaks having the best feature_type
      bestftcolnames <- intersect(
         c("fc",
            "feature_type",
            "gene_name",
            "gene_id"),
         colnames(grdt0));
      if (verbose) {
         jamba::printDebug("annotate_gr_by_genome_region(): ",
            "bestftcolnames:",
            bestftcolnames);
      }
      
      ## subset for best feature_type for each peak
      # grdt0_hasbestft <- unique(
      #    subset(grdt0,
      #       grdt0[["feature_type"]] == grdt0_bestft[grdt0[["fc"]]])[, ..bestftcolnames]);
      ## alternate syntax to avoid "..bestftcolnames"
      grdt0_hasbestft <- unique(
         subset(grdt0,
            grdt0[["feature_type"]] == grdt0_bestft[grdt0[["fc"]]]
            )[, bestftcolnames, with=FALSE]);
      
      # optionally sort only the subset rows by gene_name, gene_id
      if ("global" %in% sort_optimization) {
         sort_time <- system.time({
            grdt0_hasbestft <- jamba::mixedSortDF(grdt0_hasbestft,
               byCols=c("gene_name", "gene_id"));
         }, gcFirst=FALSE);
         if (verbose) {
            jamba::printDebug("",
               "elapsed global sort_time (sec): ", sort_time[["elapsed"]]);
         }
      }
      
      ## Comma-delimit values for each column
      i_set <- intersect(
         c("feature_type",
            "gene_name",
            "gene_id"),
         colnames(grdt0_hasbestft));
      for (i in i_set) {
         inewcolname <- paste0(i, "_winner");
         if (verbose > 1) {
            jamba::printDebug("annotate_gr_by_genome_region(): ",
               "processing inewcolname: ", inewcolname);
         }
         ikeepcolnames <- c("fc", i);
         ## avoid using "..ikeepcolnames"
         # grdt0_hasbestft_sub <- unique(grdt0_hasbestft[, ..ikeepcolnames]);
         grdt0_hasbestft_sub <- unique(
            grdt0_hasbestft[, ikeepcolnames, with=FALSE]);
         
         if (anyDuplicated(grdt0_hasbestft_sub[["fc"]])) {
            ## split only duplicate entries
            fcdupes <- is_duplicate(grdt0_hasbestft_sub[["fc"]]);
            if ("fast" %in% sort_optimization) {
               sort_time <- system.time({
                  ## avoid using "..i"
                  # grdt0_hasbestft_sub_dupe <- jamba::mixedSortDF(
                  #    grdt0_hasbestft_sub[fcdupes, c(..i, "fc"), drop=FALSE],
                  #    byCols=c(i));
                  ifc <- c(i, "fc")
                  grdt0_hasbestft_sub_dupe <- jamba::mixedSortDF(
                     grdt0_hasbestft_sub[fcdupes, ifc, with=FALSE],
                     byCols=c(i));
               }, gcFirst=FALSE);
               if (verbose) {
                  jamba::printDebug("",
                     "elapsed fast sort_time (sec): ", sort_time[["elapsed"]]);
               }
               grdt0_bestvalue_dupe <- S4Vectors::unstrsplit(
                  sep=",",
                  IRanges::CharacterList(
                     split(as.character(grdt0_hasbestft_sub_dupe[[i]]),
                        grdt0_hasbestft_sub_dupe[["fc"]])));
            } else {
               grdt0_bestvalue_dupe <- S4Vectors::unstrsplit(
                  sep=",",
                  IRanges::CharacterList(
                     split(as.character(grdt0_hasbestft_sub[[i]])[fcdupes],
                        grdt0_hasbestft_sub[["fc"]][fcdupes])));
            }
            grdt0_bestvalue_nondupe <- jamba::nameVector(
               as.character(grdt0_hasbestft_sub[[i]][!fcdupes]),
               grdt0_hasbestft_sub[["fc"]][!fcdupes]);
            grdt0_bestvalue <- c(grdt0_bestvalue_dupe,
               grdt0_bestvalue_nondupe);
         } else {
            grdt0_bestvalue <- jamba::nameVector(
               as.character(grdt0_hasbestft_sub[[i]]),
               grdt0_hasbestft_sub[["fc"]]);
         }
         imatch <- match(names(grdt0_bestvalue),
            GenomicRanges::values(gr)[[name_colname]]);
         if ("feature_type" %in% i) {
            GenomicRanges::values(gr)[, inewcolname] <- "extragenic";
         } else {
            GenomicRanges::values(gr)[, inewcolname] <- c(NA, "")[1];
         }
         if (length(imatch) > 0) {
            GenomicRanges::values(gr)[imatch, inewcolname] <- as.character(
               grdt0_bestvalue);
         }
      }
   }
   
   ##################################################
   ## gene distance analysis
   if ("nearest_gene" %in% include_type) {
      if (verbose) {
         jamba::printDebug("annotate_gr_by_genome_region(): ",
            "processing nearest_gene");
      }
      fc_gr_genedist <- GenomicRanges::distanceToNearest(
         gr,
         ignore.strand=TRUE,
         genome_regions,
         select="all");
      qh2 <- S4Vectors::queryHits(fc_gr_genedist);
      sh2 <- S4Vectors::subjectHits(fc_gr_genedist);
      keep_genecols <- intersect(
         c(gene_name_colname,
            gene_id_colname),
         colnames(GenomicRanges::values(genome_regions)));
      nearest_keep_genecols <- head(
         paste0("nearest_", c("gene_name", "gene_id")),
         length(keep_genecols));
      fc_gr_genedist_df <- data.frame(
         check.names=FALSE,
         stringsAsFactors=FALSE,
         name=GenomicRanges::values(gr[qh2])[[name_colname]],
         nearest_gene_distance=GenomicRanges::values(
            fc_gr_genedist)$distance,
         jamba::renameColumn(
            as.data.frame(GenomicRanges::values(
               genome_regions[sh2])[, keep_genecols, drop=FALSE]),
            from=keep_genecols,
            to=nearest_keep_genecols)
         );
      
      for (i in c(nearest_keep_genecols, "nearest_gene_distance")) {
         #printDebug(i);
         fc_gr_genedist_df1 <- unique(fc_gr_genedist_df[, c(name_colname, i)]);
         if (jamba::igrepHas("_distance", i)) {
            fc_gr_genedist_df1vals <- jamba::nameVector(
               fc_gr_genedist_df1[, c(i, name_colname), drop=FALSE]);
            GenomicRanges::values(gr)[, i] <- c(NA, 0)[1];
         } else {
            fc_gr_genedist_df1vals <- S4Vectors::unstrsplit(sep=",",
               IRanges::CharacterList(
                  split(fc_gr_genedist_df1[[i]],
                     fc_gr_genedist_df1[[name_colname]])
               )
            )
            GenomicRanges::values(gr)[, i] <- c(NA, "")[1];
         }
         imatch <- match(
            names(fc_gr_genedist_df1vals),
            GenomicRanges::values(gr)[[name_colname]]);
         if (length(imatch) > 0) {
            if (jamba::igrepHas("_distance", i)) {
               GenomicRanges::values(gr)[imatch, i] <- fc_gr_genedist_df1vals;
            } else {
               GenomicRanges::values(gr)[imatch, i] <- as.character(
                  fc_gr_genedist_df1vals);
            }
         }
      }
   }
   
   ##################################################
   ## Optional mask regions
   if (length(mask_regions) > 0) {
      if (verbose) {
         jamba::printDebug("annotate_gr_by_genome_region(): ",
            "processing mask_regions");
      }
      if ("character" %in% class(mask_regions)) {
         mask_regions_grl <- GenomicRanges::GRangesList(
            lapply(mask_regions, function(mask_region){
            if (file.exists(mask_region)) {
               if (verbose) {
                  jamba::printDebug("annotate_gr_by_genome_region(): ",
                     "Importing mask_region from:",
                     mask_region);
               }
               mgr <- rtracklayer::import(mask_region);
               mgr <- GenomeInfoDb::keepSeqlevels(mgr,
                  intersect(GenomeInfoDb::seqlevels(gr),
                     GenomeInfoDb::seqlevels(mgr)),
                  pruning.mode="coarse");
            }
         }));
         mask_regions_gr <- GenomicRanges::sort(mask_regions_grl@unlistData);
      }
      if ("GRangesList" %in% class(mask_regions_gr)) {
         mask_regions_gr <- mask_regions_gr@unlistData;
      }
      if ("GRanges" %in% class(mask_regions_gr)) {
         GenomicRanges::values(gr)$mask <- IRanges::overlapsAny(gr,
            mask_regions_gr);
      }
   }

   return(gr)   
}


#' Flatten genome_regions
#' 
#' Flatten genome_regions
#' 
#' This function is a simple wrapper method to take `genome_regions`
#' and produce a flattened version with just the feature_type winners.
#' It calls `annotate_gr_from_genome_regions()` in a somewhat
#' streamlined way.
#' 
#' The chromosome lengths are either used from the `genome_regions` object,
#' or by using argument `genome` to query UCSC data via
#' `GenomeInfoDb::fetchExtendedChromInfoFromUCSC()`. These chromosome lengths
#' are used to create chromosome-length `GRanges`, so that any region not
#' described by `genome_regions` will be considered `extragenic`. In this
#' way, every region in the genome will be represented.
#' 
#' @family slicejam genome regions
#' 
#' @param genome_regions `GRanges` object as produced by
#'    `genomic_regions_from_gtf()`
#' @param genome `character` string indicating the genome, recognized
#'    by `GenomeInfoDb::fetchExtendedChromInfoFromUCSC()` in order to
#'    obtain chromosome lengths. This value is only used when
#'    `seqlengths(genome_regions)` contains `NA` values.
#' @param name_colname,gene_name_colname,feature_type_colname,gene_id_colname
#'    `character` strings indicating colnames in `colnames(values(genome_regions))`
#'    and passed to `annotate_gr_by_genome_region()`. Note that only
#'    `feature_type_colname` and `gene_name_colname` are included in output.
#' @param feature_grep_order `character` vector passed to
#'    `annotate_gr_by_genome_region()` to define priority of feature type
#'    values. The value is ultimately passed to `jamba::provigrep()`.
#' @param mask_regions passed to `annotate_gr_by_genome_region()` and
#'    is currently ignored for this function.
#' @param canonical_only `logical` indicating whether to filter `genome_regions`
#'    to include only canonical chromosome names, for example starting with
#'    `"chr"` and including no underscore `"_"`, period `"."` or hyphen/dash
#'    `"-"` characters.
#' @param sort_optimization `character` string indicating the method of sorting:
#'    * `"fast"` sorts `gene_name` independently from `gene_id`
#'    * `"global"` sorts `gene_name` then `gene_id`, and maintains this order
#' @param verbose `logical` indicating whether to print verbose output
#' @param ... additional arguments are ignored.
#' 
#' @export
flatten_genome_regions <- function
(genome_regions,
 genome="hg19",
 name_colname="name",
 gene_name_colname="gene_name",
 feature_type_colname="feature_type",
 gene_id_colname="gene_id",
 feature_grep_order = c("promoter",
    "tts",
    "exon",
    "intron",
    "extragenic",
    "intergenic",
    "."),
 mask_regions=NULL,
 canonical_only=TRUE,
 sort_optimization="fast",
 verbose=FALSE,
 ...)
{
   # if (!suppressPackageStartupMessages(require(GenomicRanges))) {
   #    stop("GenomicRanges package is required.");
   # }
   # if (!suppressPackageStartupMessages(require(GenomeInfoDb))) {
   #    stop("GenomeInfoDb package is required.");
   # }
   
   # check non-canonical chromosomes
   is_noncanonical_gr <- grepl("^[^c][^h][^r]|[-._]",
      ignore.case=TRUE,
      GenomeInfoDb::seqlevels(genome_regions));
   if (canonical_only && any(is_noncanonical_gr)) {
      # enforce seqlevels
      canonical_seqlevels <- GenomeInfoDb::seqlevels(
         genome_regions)[!is_noncanonical_gr];
      GenomeInfoDb::seqlevels(genome_regions,
         pruning.mode="coarse") <- canonical_seqlevels;
   }
   
   # get chromosome sizes
   if (all(is.na(GenomeInfoDb::seqlengths(genome_regions)))) {
      chrom_info <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(genome);
      if (FALSE && canonical_only) {
         is_noncanonical_ci <- grepl("^[^c][^h][^r]|[-._]",
            ignore.case=TRUE,
            chrom_info$UCSC_seqlevel);
         chrom_info <- subset(chrom_info, !is_noncanonical_ci);
      }
      pmax_na <- function(...){pmax(..., na.rm=TRUE)}
      seq_match <- Reduce("pmax_na",
         list(
            match(
               GenomeInfoDb::seqlevels(genome_regions),
               chrom_info$UCSC_seqlevel),
            match(
               GenomeInfoDb::seqlevels(genome_regions),
               chrom_info$GenBankAccn)));
      seqlengths(genome_regions) <- chrom_info$UCSC_seqlength[seq_match];
      genome(genome_regions) <- genome;
      isCircular(genome_regions) <- chrom_info$circular[seq_match];
   }

   # make chromosome-length features
   chrom_gr <- GenomicRanges::GRanges(
      seqnames=names(GenomeInfoDb::seqlengths(genome_regions)),
      ranges=IRanges::IRanges(start=1,
         end=GenomeInfoDb::seqlengths(genome_regions)),
      strand="*");
   extragenic_gr <- GenomicRanges::setdiff(chrom_gr,
      genome_regions,
      ignore.strand=TRUE);
   
   # flatten input regions
   genome_regions_flat <- GenomicRanges::disjoin(genome_regions,
      ignore.strand=TRUE);
   if (length(extragenic_gr) > 0) {
      genome_regions_flat <- sort(c(genome_regions_flat,
         extragenic_gr))
   }
   GenomicRanges::values(genome_regions_flat)$name <- paste0("region",
      jamba::padInteger(seq_along(genome_regions_flat)));
   
   # annotate using same method used for peaks
   genome_regions_ann <- annotate_gr_by_genome_region(
      gr=genome_regions_flat,
      genome_regions=genome_regions[, c(gene_name_colname,
         feature_type_colname)],
      name_colname=name_colname,
      gene_name_colname=gene_name_colname,
      gene_id_colname=gene_id_colname,
      feature_type_colname=feature_type_colname,
      mask_regions=mask_regions,
      feature_grep_order=feature_grep_order,
      include_type="winner",
      sort_optimization=sort_optimization,
      verbose=verbose,
      ...);
   
   # adjust feature_type_winner
   GenomicRanges::values(genome_regions_ann)$feature_type_winner <- factor(
      GenomicRanges::values(genome_regions_ann)$feature_type_winner,
      levels=jamba::provigrep(
         c("promoter", "tts", "exon", "intron", "extra", "."),
         unique(GenomicRanges::values(genome_regions_ann)$feature_type_winner))
      );
   
   return(genome_regions_ann);
}

#' Convert promoter region to a single TSS position
#' 
#' Convert promoter region to a single TSS position
#' 
#' This function is called by `slicejam_analysis.Rmd`. It takes the promoter
#' region defined around a transcript start site (TSS) by `upstream` and
#' `downstream` distances from the TSS, and instead returns the original
#' TSS position.
#' 
#' ### Note
#' 
#' This function detects regions which are smaller than
#' `upstream + downstream`, which would cause an error during adjustment.
#' In this case, features at the start of the chromosome with start `0` or `1`
#' will adjust the `GenomicRanges::end()` position, then will define the
#' start position to have `width=1`.
#' 
#' In that case this function should trust the point farthest from the end
#' of the chromosome, since that distance would not have been truncated by
#' the chromosome size.
#' 
#' This step makes the assumption of the expected
#' starting region width, which may not be correct and may need to be
#' an optional argument. For example `upstream=2000,downstream=200` may
#' imply a region with `width=2200`, or `width=2201`, depending upon
#' the source data.
#' 
#' Bonus points: calculate `width(gr) - upstream - downstream` and take
#' the most frequently occuring value. If this value is `0` or `1` then
#' proceed to check for peaks that are smaller.
#' 
#' 
#' @return `GRanges` object after adjusting the `upstream` and `downstream`
#'    positions, relative to the `GenomicRanges::strand()` of each feature.
#' 
#' @family slicejam genome regions
#' 
#' @param gr `GRanges` object representing TSS regions defined using the
#'    `upstream` and `downstream` argument values.
#' @param upstream,downstream `integer` number of bases upstream, and downstream
#'    the stranded gene TSS position, used to define
#'    a "promoter region" around that site.
#' @param minimum_width `integer` width below which a feature will be
#'    adjusted only using the start or end of the range, and the output
#'    will automatically have `width=1`.
#' @param positive `character` vector indicating `GenomicRanges::strand()`
#'    values processed as positive strands.
#' 
#' @export
promoter_to_tss <- function
(gr,
 upstream=2000,
 downstream=200,
 minimum_width=upstream + downstream,
 positive=c("+", "*"))
{
   ## Purpose is to convert coordinates of promoters back to original TSS
   to_process <- (GenomicRanges::width(gr) >= (minimum_width));
   
   
   # Process entries with at least the minimum width, on "+" or "-" strand
   gr2pos <- (as.character(GenomicRanges::strand(gr)) %in% positive);
   newstart <- ifelse(gr2pos[to_process],
      GenomicRanges::start(gr[to_process]) + upstream,
      GenomicRanges::start(gr[to_process]) + downstream - 1);
   newend <- ifelse(gr2pos[to_process],
      GenomicRanges::end(gr[to_process]) - downstream + 1,
      GenomicRanges::end(gr[to_process]) - upstream);
   GenomicRanges::ranges(gr[to_process]) <- IRanges::IRanges(
      start=newstart,
      end=newend);
   
   # Handle entries too small for processing above
   if (any(!to_process)) {
      front_clip <- GenomicRanges::start(gr[!to_process]) <= 1;
      # A front_clip happens when one end of the feature is at
      # the front of a chromosome, suggesting it was clipped here.
      # Otherwise assume it was clipped at the far end of the chromosome.
      # For front_clip trust the end(),
      # otherwise trust the start().
      newstart_small <- ifelse(front_clip,
         ifelse(gr2pos[!to_process],
            GenomicRanges::end(gr[!to_process]) - downstream,
            GenomicRanges::end(gr[!to_process]) - upstream),
         ifelse(gr2pos[!to_process],
            GenomicRanges::start(gr[!to_process]) + upstream,
            GenomicRanges::start(gr[!to_process]) + downstream))
      newend_small <- ifelse(front_clip,
         ifelse(gr2pos[!to_process],
            (GenomicRanges::end(gr[!to_process]) - downstream) + 1,
            (GenomicRanges::end(gr[!to_process]) - upstream) + 1),
         ifelse(gr2pos[!to_process],
            GenomicRanges::start(gr[!to_process]) + upstream + 1,
            GenomicRanges::start(gr[!to_process]) + downstream + 1));
      GenomicRanges::ranges(gr[!to_process]) <- IRanges::IRanges(
         start=newstart_small,
         end=newend_small);
   }
   return(gr);
}

#' Trim multiple gene descriptions to max number of entries per field
#' 
#' Trim multiple gene descriptions to max number of entries per field
#' 
#' This function is internally called by `slicejam_analysis.Rmd` and is
#' not intended for much wider use.
#' 
#' @param df `data.frame`
#' @param max_genes `integer` maximum number of delimited entries to keep.
#' @param symbol_colname `character` indicating the column to test for
#'    multiple entries.
#' @param gene_colnames `character` indicating the columns to update
#'    for rows where multiple entries are detected in `symbol_colname`.
#' @param sep `character` string of the expected delimiter between multiple
#'    values.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @export
trim_multigenedesc <- function
(df,
 max_genes=4,
 symbol_colname="SYMBOL",
 gene_colnames=c("SYMBOL", "GENENAME"),
 sep=";",
 verbose=FALSE,
 ...)
{
   # build grap patterns
   if (!symbol_colname %in% colnames(df)) {
      stop("symbol_colname must be present in colnames(df)")
   }
   if (!any(gene_colnames %in% colnames(df))) {
      stop("gene_colnames must be present in colnames(df)")
   }
   gene_colnames <- intersect(gene_colnames,
      colnames(df));
   if (!length(sep) == 1) {
      stop("sep must have length=1")
   }
   multigrep <- paste0("^(",
      paste(
         rep(
            paste0("[^", sep, "]+"),
            max_genes),
         collapse=";"),
      ");.+$");
   
   # detect rows to process
   multirows <- grepl(multigrep, df[[symbol_colname]]);
   if (any(multirows)) {
      if (verbose) {
         jamba::printDebug(c("   Trimming ",
            jamba::formatInt(sum(multirows)),
            " multi-gene rows."),
            sep="");
      }
      for (icol in gene_colnames) {
         df[[icol]] <- gsub(multigrep,
            "\\1",
            df[[icol]]);
      }
   }
   return(df);
}

