
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
#' exon, or intron, is referred to as `"intergenic"`, although intergenic
#' regions are not defined by this function directly.
#' 
#' In more detail: The arguments `geneAttrNames`, `txAttrNames`, and
#' `geneFeatureType` are used by `splicejam::makeTx2geneFromGtf()` to
#' create a `data.frame` which relates `transcript_id` to each `gene_name`
#' by default. If the GTF file uses different attributes, adjust
#' `geneAttrNames` and `txAttrNames` to match the GTF file.
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
      refgene_txdb <- GenomicFeatures::makeTxDbFromGFF(gtf);
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
   # rename colnames
   GenomicRanges::values(txByGene@unlistData) <- jamba::renameColumn(
      GenomicRanges::values(txByGene@unlistData),
      from=c("tx_id", "tx_name"),
      to=c("internal_tx_id", tx_colname));
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
   values(ttsByTx) <- values(txByGene@unlistData);
   # tts reduce() per gene
   ttsByGene <- GenomicRanges::split(ttsByTx,
      GenomicRanges::values(ttsByTx)[[gene_colname]]);
   ttsByGeneRed <- GenomicRanges::reduce(ttsByGene);
   GenomicRanges::values(ttsByGeneRed@unlistData)[[gene_colname]] <- rep(
      names(ttsByGeneRed),
      lengths(ttsByGeneRed));
   gene_match <- match(values(ttsByGeneRed@unlistData)[[gene_colname]],
      values(ttsByGene@unlistData)[[gene_colname]])
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
   values(promoters_gr) <- values(txByGene@unlistData);
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
      gene_match <- match(values(tssByGeneRed@unlistData)[[gene_colname]],
         values(tssByGene@unlistData)[[gene_colname]])
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
               provigrep(c("gene.*name", "."),
                  c(geneAttrNames, txAttrNames)), 
               1),
            "feature_type");
         gr_id <- paste0(seqnames(genome_regions), ":", start(genome_regions), "-", end(genome_regions));
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
 feature_grep_order=c("promoter", "tts", "exon", "intron", "intergenic", "."),
 include_type=c("overlap", "winner", "nearest_gene"),
 verbose=FALSE,
 ...)
{
   # validate input
   include_type <- match.arg(include_type,
      several.ok=TRUE);
   
   name_colname <- intersect(name_colname,
      colnames(values(gr)));
   if (length(name_column) == 0) {
      name_column <- "name";
      values(gr)[[name_column]] <- paste0("gr",
         jamba::padInteger(seq_along(gr)));
   }
   
   gene_name_colname <- head(intersect(gene_name_colname,
      colnames(values(genome_regions))), 1);
   feature_type_colname <- head(intersect(feature_type_colname,
      colnames(values(genome_regions))), 1);
   if (length(gene_name_colname) == 0 || length(feature_type_colname) == 0) {
      stop("gene_name_colname and feature_type_colname must be present in colnames(values(genome_regions))");
   }
   gene_id_colname <- intersect(gene_id_colname,
      colnames(values(genome_regions)));
   
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
      grdt <- data.table::data.table(
         fc=GenomicRanges::values(gr[S4Vectors::queryHits(fco)])[[name_colname]],
         gene_name=GenomicRanges::values(genome_regions_exp[S4Vectors::subjectHits(fco)])[[gene_name_colname]],
         #gene_id=values(genome_regions_exp[subjectHits(fco)])$gene_id,
         feature_type=factor(
            GenomicRanges::values(genome_regions_exp[S4Vectors::subjectHits(fco)])[[feature_type_colname]],
            levels=provigrep(feature_grep_order,
               unique(GenomicRanges::values(genome_regions_exp[S4Vectors::subjectHits(fco)])[[feature_type_colname]]))),
         gene_feature_type=GenomicRanges::values(genome_regions_exp[S4Vectors::subjectHits(fco)])$gene_feature_type);
      if (FALSE) {
         grdt <- jamba::renameColumn(grdt,
            from=c("gene_name",
               "feature_type"),
            to=c(gene_name_colname,
               feature_type_colname));
      }

      if (length(gene_id_colname) > 0) {
         grdt$gene_id <- GenomicRanges::values(genome_regions_exp[S4Vectors::subjectHits(fco)])[[gene_id_colname]];
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
            idt <- grdt[, ..icols];
            ## expand comma-delimited entries
            grdt0ft <- unique(jamba::mixedSortDF(idt, byCols=1:2)[,c(1,2)]);
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
            if (feature_type_colname == i) {
               GenomicRanges::values(gr)[,i] <- unname(c("intergenic", grd_vals[[i]][1])[1]);
            } else {
               GenomicRanges::values(gr)[,i] <- unname(c(NA, grd_vals[[i]][1])[1]);
            }
            GenomicRanges::values(gr)[imatch,i] <- grd_vals[[i]];
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
      grdt0 <- jamba::mixedSortDF(grdt,
         byCols=c("fc",
            "feature_type",
            "gene_name",
            "gene_id"));
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
      
      grdt0_hasbestft <- unique(
         subset(grdt0,
            grdt0[["feature_type"]] == grdt0_bestft[grdt0[["fc"]]])[, ..bestftcolnames]);
      
      ## Comma-delimit values for each column
      ## peak_089268
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
         grdt0_hasbestft_sub <- unique(grdt0_hasbestft[, ..ikeepcolnames]);
         if (anyDuplicated(grdt0_hasbestft_sub[["fc"]])) {
            ## split only duplicate entries
            fcdupes <- is_duplicate(grdt0_hasbestft_sub[["fc"]]);
            grdt0_bestvalue_dupe <- S4Vectors::unstrsplit(
               sep=",",
               IRanges::CharacterList(
                  split(grdt0_hasbestft_sub[[i]][fcdupes],
                     grdt0_hasbestft_sub[["fc"]][fcdupes])));
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
         if (feature_type_colname %in% i) {
            GenomicRanges::values(gr)[,inewcolname] <- "intergenic";
         } else {
            GenomicRanges::values(gr)[,inewcolname] <- c(NA, "")[1];
         }
         if (length(imatch) > 0) {
            GenomicRanges::values(gr)[imatch,inewcolname] <- as.character(grdt0_bestvalue);
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
         name=GenomicRanges::values(gr[S4Vectors::queryHits(fc_gr_genedist)])[[name_colname]],
         nearest_gene_distance=GenomicRanges::values(fc_gr_genedist)$distance,
         jamba::renameColumn(
            as.data.frame(GenomicRanges::values(genome_regions[S4Vectors::subjectHits(fc_gr_genedist)])[,keep_genecols]),
            from=keep_genecols,
            to=nearest_keep_genecols)
         );
      
      for (i in c(nearest_keep_genecols, "nearest_gene_distance")) {
         #printDebug(i);
         fc_gr_genedist_df1 <- unique(fc_gr_genedist_df[,c(name_colname, i)]);
         if (jamba::igrepHas("_distance", i)) {
            fc_gr_genedist_df1vals <- jamba::nameVector(fc_gr_genedist_df1[,c(i, name_colname)]);
            GenomicRanges::values(gr)[,i] <- c(NA, 0)[1];
         } else {
            fc_gr_genedist_df1vals <- S4Vectors::unstrsplit(sep=",",
               IRanges::CharacterList(
                  split(fc_gr_genedist_df1[[i]],
                     fc_gr_genedist_df1[[name_colname]])
               )
            )
            GenomicRanges::values(gr)[,i] <- c(NA, "")[1];
         }
         imatch <- match(
            names(fc_gr_genedist_df1vals),
            GenomicRanges::values(gr)[[name_colname]]);
         if (length(imatch) > 0) {
            if (jamba::igrepHas("_distance", i)) {
               GenomicRanges::values(gr)[imatch,i] <- fc_gr_genedist_df1vals;
            } else {
               GenomicRanges::values(gr)[imatch,i] <- as.character(fc_gr_genedist_df1vals);
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
         mask_regions_grl <- GenomicRanges::GRangesList(lapply(mask_regions, function(mask_region){
            if (file.exists(mask_region)) {
               if (verbose) {
                  jamba::printDebug("annotate_gr_by_genome_region(): ",
                     "Importing mask_region from:",
                     mask_region);
               }
               mgr <- rtracklayer::import(mask_region);
               mgr <- GenomeInfoDb::keepSeqlevels(mgr,
                  intersect(GenomeInfoDb::seqlevels(gr),
                     GenomeInfoDb::seqlevels(mgr)));
            }
         }));
         mask_regions_gr <- GenomicRanges::sort(mask_regions_grl@unlistData);
      }
      if ("GRangesList" %in% class(mask_regions_gr)) {
         mask_regions_gr <- mask_regions_gr@unlistData;
      }
      if ("GRanges" %in% class(mask_regions_gr)) {
         GenomicRanges::values(gr)$mask <- IRanges::overlapsAny(gr, mask_regions_gr);
      }
   }

   return(gr)   
}
