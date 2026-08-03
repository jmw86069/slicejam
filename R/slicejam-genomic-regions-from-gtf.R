
#' Genomic Regions from GTF and BED
#' 
#' Genomic Regions from GTF and BED
#' 
#' This function represents the 'Rscript' functionality
#' of 'genomic_regions_from_gtf', by wrapping:
#' 
#' * `genomic_regions_from_gtf()` which creates the
#' genomic regions from GTF, using detected genes and
#' transcripts when provided, and using custom promoter
#' and TTS region sizes.
#' * `annotate_gr_by_genome_region()` which uses the
#' BED file data (GRanges) and annotates them with
#' the genomic regions.
#' 
#' Note the Rscript will prepare genome regions when
#' there is no BED file. If you have no BED file to use,
#' use `genomic_regions_from_gtf()` which does not use
#' the BED file.
#' 
#' @returns `list` with supporting data objects.
#' The key object is 'bed_gr' containing the BED data
#' annotated with genome regions.
#' * 'genome_regions': `GRanges` with genomic regions
#' as defined by the other arguments.
#' * 'mask': `GRanges` only when provided as a BED file
#' or `GRanges` object.
#' * 'bed_gr': `GRanges` of 'bed' data annotated using
#' the genomic regions.
#' * 'outfile': the output file when 'save_outdata=TRUE'.
#' * 'outfile_rdata': the output RData file when
#' 'save_rdata=TRUE'.
#' 
#' @family slicejam genome regions
#' 
#' @param gtf `character` path to GTF file
#' @param bed `character` path to BED file, or `GRanges`.
#' @param upstream_promoter,downstream_promoter `integer` distance
#'    upstream or downstream the gene transcript start site (TSS),
#'    with default 2000, 200, respectively.
#' @param upstream_tts,downstream_TTS `integer` distance
#'    upstream or downstream the gene transcript termination site
#'    (TTS), with default 1000, 1000, respectively.
#' @param detectedTx `character` path to file with detected
#'    transcripts, which match 'transcript_id' in the GTF file.
#' @param detectedGenes `character` path to file with detected
#'    genes, which match 'gene_name' in the GTF file.
#' @param mask `character` path to BED file with mask regions, or
#'    `GRanges` object.
#' @param save_outfile `logical`, default FALSE, whether to save
#'    the results to an annotated BED file.
#' @param outfile `character` output file, default NULL uses the BED
#'    filename to create a filename with `'.genome_regions.'`, the
#'    current date, and ending `'.txt'`.
#' @param save_rdata `logical` default FALSE, whether to save data
#'    as an RData file. It uses the `'outfile'` and replaces 'txt'
#'    with 'RData'.
#' @param rdata_file `character` default NULL, to specify a specific
#'    RData output file.
#' @param save_txdb `logical` default TRUE, whether to save the
#'    `TxDb` created from the GTF file. This step saves time when
#'    using the same GTF in future.
#' @param force_refresh `logical` default FALSE, forces the GTF
#'    to be reprocessed, therefore not re-using any `TxDb` data.
#' @param max_threads `integer` threads used by 'data.table',
#'    default 4.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @export
genomic_regions_from_gtf_bed <- function
(gtf,
 bed=NULL,
 upstream_promoter=2000,
 downstream_promoter=200,
 upstream_tts=1000,
 downstream_tts=1000,
 detectedTx=NULL,
 detectedGenes=NULL,
 mask=NULL,
 save_outfile=FALSE,
 outfile=NULL,
 save_rdata=FALSE,
 rdata_file=NULL,
 save_txdb=TRUE,
 force_refresh=FALSE,
 max_threads=4,
 verbose=FALSE,
 ...)
{
   #
   retvals <- list();

   # check detectedGenes as file
   if (length(detectedGenes) > 0 && all(file.exists(detectedGenes))) {
      detectedGenes <- unname(unique(unlist(
         lapply(detectedGenes, readLines))));
   }

   # check detectedTx as file
   if (length(detectedTx) > 0 && all(file.exists(detectedTx))) {
      detectedTx <- unname(unique(unlist(
         lapply(detectedTx, readLines))));
   }

   # Set data.table threads
   data.table::setDTthreads(threads=max_threads)

   if (verbose) {
      jamba::printDebug("Confirming genomic regions data.");
   }
   st1 <- system.time({
      genome_regions <- genomic_regions_from_gtf(gtf=gtf,
         upstream_promoter=upstream_promoter,
         downstream_promoter=downstream_promoter,
         upstream_tts=upstream_tts,
         downstream_tts=downstream_tts,
         detectedTx=detectedTx,
         detectedGenes=detectedGenes,
         save_rdata=save_rdata,
         rdata_file=rdata_file,
         save_txdb=save_txdb,
         force_refresh=force_refresh,
         verbose=verbose);
   })
   if (verbose) {
      jamba::printDebug("Genomic regions time elapsed: ",
         unlist(st1["elapsed"]));
   }
   retvals$genome_regions <- genome_regions;

   # Load Mask BED
   if (length(mask) > 0) {
      if (inherits(mask, "GRanges")){
         # use mask as-is
      } else if (file.exists(mask)){
         mask <- rtracklayer::import(mask);
      } else {
         stop("mask must be a BED file, or GRanges.")
      }
      # sort seqlevels in proper chromosome order
      mask <- GenomeInfoDb::keepSeqlevels(mask,
         jamba::mixedSort(GenomeInfoDb::seqlevels(mask)));
      retvals$mask <- mask;
   }
   
   # Load BED data
   if (length(bed) > 0) {
      if (file.exists(bed)) {
         #bed <- "/path_to/project/BRG1KD_Strengthened_peaks_q005_FC1_5.bed";
         if (verbose) {
            jamba::printDebug("Importing bed file.");
         }
         fc_gr <- rtracklayer::import(bed);
      } else if (inherits(bed, "GRanges")) {
         fc_gr <- bed;
      }

      # make sure "name" column is present
      if (!"name" %in% colnames(GenomicRanges::values(fc_gr))) {
         GenomicRanges::values(fc_gr)[,"name"] <- paste0("bed_",
            jamba::padInteger(seq_along(fc_gr)));
      }

      # sort seqlevels in proper chromosome order
      fc_gr <- GenomeInfoDb::keepSeqlevels(fc_gr,
         jamba::mixedSort(GenomeInfoDb::seqlevels(fc_gr)));
      
      if (anyDuplicated(GenomicRanges::values(fc_gr)$name)) {
         if (verbose) {
            jamba::printDebug("Enforcing new unique peak names.");
         }
         GenomicRanges::values(fc_gr)$bed_name <- GenomicRanges::values(fc_gr)$name;
         GenomicRanges::values(fc_gr)$name <- paste0("peak_",
            jamba::padInteger(seq_along(fc_gr)));
      }
      
      ## annotate using genome_regions
      if (verbose) {
         jamba::printDebug("Annotating bed with genome regions.");
      }
      st2 <- system.time({bed_gr <- annotate_gr_by_genome_region(
         gr=fc_gr,
         genome_regions=genome_regions,
         mask_regions=mask,
         verbose=verbose > 1)});
      if (verbose) {
         jamba::printDebug("Annotation time elapsed: ",
         unlist(st2["elapsed"]));
      }
      retvals$bed_gr <- bed_gr;
      
      if (length(outfile) == 0 || nchar(outfile) == 0) {
         outfile <- gsub("[.](bed[0-9]*|broadPeak|narrowPeak|peak)$",
            ignore.case=TRUE,
            "",
            bed);
         outfile <- paste0(outfile,
            ".genome_regions.",
            jamba::getDate(),
            ".txt");
      }
      if (isTRUE(save_outfile)) {
         retvals$outfile <- outfile;
         data.table::fwrite(
            x=as.data.frame(bed_gr),
            file=outfile,
            sep="\t");
         if (verbose) {
            jamba::printDebug("Output saved: ",
               outfile);
         }      
      }
      if (isTRUE(save_rdata)) {
         outfile_rdata <- gsub("[.]txt$", ".RData", outfile);
         retvals$outfile_rdata <- outfile_rdata;
         save(list=c("bed_gr"), file=outfile_rdata);
         if (verbose) {
            jamba::printDebug("RData saved: ",
               outfile_rdata);
         }      
      }
   }
   invisible(retvals);
}
