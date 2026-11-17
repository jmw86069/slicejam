
#' Get stranded distance relative to the subject genome ranges.
#' 
#' @returns `numeric` signed distance of `gr1` features relative
#'    to `gr2` features and the `gr2` strandedness.
#' 
#' @param d `Hits` output from `GenomicRanges::distanceToNearest()`
#'    which must include 'distance' in the results.
#' @param gr1,gr2 `GenomicRanges` with query and subject, respectively.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @family slicejam genome regions
#' 
#' @keywords internal
#' @noRd
get_stranded_distance <- function(
   d,
   gr1,
   gr2,
   verbose=FALSE,
   ...
) {
   #
   if (!"distance" %in% names(GenomicRanges::values(d))) {
      stop("The 'd' object must have values with 'distance'.")
   }
   # Logic:
   # - if distance == 0, no sign is necessary
   # - if distance > 0, either the start or end of either feature should be sufficient
   dist_values <- GenomicRanges::values(d)$distance;
   k1 <- S4Vectors::queryHits(d);
   k2 <- S4Vectors::subjectHits(d);
   start1 <- GenomicRanges::start(gr1[k1]);
   start2 <- GenomicRanges::start(gr2[k2]);
   strand2 <- GenomicRanges::strand(gr2[k2]);
   diff12 <- (start1 - start2) * sign(dist_values)
   use_sign <- ifelse(
      dist_values == 0,
      0,
      ifelse(
         as.character(strand2) %in% c("+", "*"),
         ifelse(start1 < start2, -1, 1),
         ifelse(start1 < start2, 1, -1)
      )
   )
   new_dist_values <- dist_values * use_sign;
   if (verbose) {
      print(data.frame(
         dist_values,
         start1,
         start2,
         strand2,
         diff12,
         use_sign,
         new_dist_values));# debug
   }
   new_dist_values;
}

#' Find the nearest genomic range neighbor, with stranded distance
#' 
#' Find the nearest genomic range neighbor, with stranded distance
#' relative to the subject.
#' 
#' This function is a lightweight extension to method
#' `GenomicRanges::distanceToNearest()` which also calls
#' `get_stranded_distance()`, then populates the stranded distance
#' into the 'distance' column of the `Hits` object.
#' The stranded distance is "signed":
#' * `'0'`: there is no distance, the genomic range overlaps directly.
#' * `'+'`: the query `x` is downstream of the `subject`, relative to
#' the 'strand' of the genomic range in the `subject`.
#'    * For `subject` on the '+' strand, the query has higher coordinates.
#'    * For `subject` on the '-' strand, the query has lower coordinates.
#' 
#' * `'-'`: the query `x` is upnstream of the `subject`, relative to
#' the 'strand' of the genomic range in the `subject`.
#'    * For `subject` on the '+' strand, the query has lower coordinates.
#'    * For `subject` on the '-' strand, the query has higher coordinates.
#' 
#' @param x,subject `GenomicRanges` for the query and subject, respectively.
#' @param ignore.strand `logical` default TRUE, whether to match strandedness
#'    of both features. TRUE is recommended in order to report the stranded
#'    distance of "nearest match on either strand", which reports the
#'    distance relative to the subject.
#' @param ... additional arguments are passed to
#'    `GenomicRanges::distanceToNearest()`.
#' 
#' @returns `Hits` object with 'distance' containing stranded values:
#'    negative values are upstream, positive values are downstream,
#'    relative to the strandedness of the `subject`.
#' 
#' @family slicejam genome regions
#' 
#' @examples
#' genome_regions1 <- GenomicRanges::GRanges(
#'    seqnames=c("chr1"),
#'    ranges=IRanges::IRanges(
#'       start=c(1501, 3001, 3901, 5501, 6001, 6501, 6001),
#'         end=c(3500, 3900, 5500, 6000, 6500, 7000, 8000),
#'    ),
#'    strand=c("+"),
#'    gene_name="TESTGENE1",
#'    gene_id="ENSG0000001.1",
#'    feature_type=c("Promoters (-1500,+500)", "exon", "intron",
#'       "exon", "intron", "exon", "TTS (-1000,+1000)"
#'    )
#' )
#' genome_regions2 <- genome_regions1;
#' end(genome_regions2) <- end(genome_regions1) + 35000;
#' start(genome_regions2) <- start(genome_regions1) + 35000;
#' values(genome_regions2)$gene_name <- "TESTGENE2";
#' values(genome_regions2)$gene_id <- "ENSG0000002.1"
#' 
#' start3 <- 28001 - end(genome_regions1);
#' end3 <- 28001 - start(genome_regions1);
#' genome_regions3 <- GenomicRanges::GRanges(
#'    seqnames=c("chr1"),
#'    ranges=IRanges::IRanges(start=rev(start3), end=rev(end3)),
#'    strand="-",
#'    gene_name="TESTGENE3",
#'    gene_id="ENSG0000003.1",
#'    feature_type=rev(values(genome_regions1)$feature_type)
#' )
#' 
#' genome_regions <- c(genome_regions1, genome_regions2, genome_regions3)
#' 
#' peaks <- GenomicRanges::GRanges(
#'    seqnames=c("chr1"),
#'    ranges=IRanges::IRanges(
#'       start=c(1201, 1801, 3841, 4501, 10001, 19001, 30001, 36001),
#'         end=c(1300, 1900, 3940, 4600, 10100, 19100, 30100, 36100),
#'    ),
#'    strand=c("+"),
#'    gene_name=paste0("peak_00", 1:8),
#'    gene_id=paste0("peaks"),
#'    feature_type="peak"
#' )
#' genome_regions_sub <- subset(genome_regions, !grepl("Promoter|TTS", feature_type))
#' dtn <- distanceToNearestStranded(peaks, genome_regions_sub)
#' values(peaks)$nearest_gene <- genome_regions_sub$gene_name[subjectHits(dtn)]
#' values(peaks)$nearest_gene_strand <- as.character(strand(genome_regions_sub)[subjectHits(dtn)])
#' values(peaks)$nearest_range <- paste0(start(genome_regions_sub), "-", end(genome_regions_sub))[subjectHits(dtn)]
#' values(peaks)$nearest_feature_type <- genome_regions_sub$feature_type[subjectHits(dtn)]
#' values(peaks)$nearest_distance <- values(dtn)$distance
#' peaks
#' 
#' @export
distanceToNearestStranded <- function(
   x,
   subject=NULL,
   ignore.strand=TRUE,
   ...
) {
   #
   dtn <- GenomicRanges::distanceToNearest(
      x=x,
      subject=subject,
      ignore.strand=ignore.strand,
      ...)
   if (length(subject) == 0) {
      subject <- x;
   }
   distance_stranded <- get_stranded_distance(
      d=dtn,
      gr1=x,
      gr2=subject
   );
   GenomicRanges::values(dtn)$distance <- distance_stranded;
   return(dtn)
}
