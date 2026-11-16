
#' Convert promoter or TTS region to a single TSS or TTS position
#' 
#' Convert promoter region to a single TSS position, or
#' convert TTS region to a single TTS position.
#' For example, reverse the process of applying upstream
#' and downstream width to define a 'promoter' or 'TTS' region,
#' based upon name format 'Promoters (-2000, +200)' or
#' 'TTS (-1000, +1000)'.
#' 
#' This function is called by `slicejam_analysis.Rmd`, and
#' by functions `genomic_regions_by_gtf()`,
#' `genomic_regions_by_gtf_bed()` to determine nearest gene TSS.
#' 
#' It takes the promoter region defined around a
#' transcript start site (TSS) by `upstream` and
#' `downstream` distances from the TSS, and instead returns
#' the original TSS position.
#' 
#' ### Note
#' 
#' This function will also adjust for positions near the
#' chromosome end, where for example part of the region
#' 'Promoters (-2000, +-200)' is beyond the chromosome.
#' It handles this situation by using whichever side of
#' the feature is not beyond the chromosome.
#' 
#' This step makes the assumption of the expected
#' starting region width, which may not be correct depending
#' how the original operation was performed (if outside of
#' slicejam).
#' For example 'upstream=2000' and 'downstream=200' should
#' define a region with `width=2200`, however it could also
#' be `width=2201`. Use the appropriate value with
#' argument 'minimum_width'.
#' 
#' ## Bonus points
#' 
#' Calculate `width(gr) - upstream - downstream` and take
#' the most frequently occuring value. If this value is not
#' consistently `0` or `1` then it should be reviewed. A mixture
#' of `0` and `1` indicates the TSS may not have been defined
#' with consistent rules.
#' 
#' @returns `GRanges` object after adjusting the `upstream` and
#'    `downstream` ranges, relative to the
#'    `GenomicRanges::strand()` of each feature.
#' 
#' @family slicejam genome regions
#' 
#' @param gr `GRanges` object representing TSS regions defined using the
#'    `upstream` and `downstream` argument values.
#' @param use_feature_type `logical` whether to use the 'feature_type'
#'    column to populate the upstream and downstream distances.
#' @param upstream,downstream `integer` number of bases upstream, and
#'    downstream the stranded gene TSS position, used to define
#'    a "promoter region" around that site.
#'    Default is 2000, 200, respectively, but is superceded when input
#'    `GRanges` has 'feature_type' column with proper numeric ranges,
#'    for example 'Promoter (-2000,+200)' or 'TSS(-1000,+1000)'.
#' @param minimum_width `integer`,
#'    width below which a feature is considered to be truncated,
#'    to be adjusted only using the start or end of the range.
#'    Default is the sum of 'upstream' and 'downstream'.
#' @param tts_upstream,tts_downstream `numeric` default 1000, applied
#'    to TSS upstream and downstream regions.
#'    These values are superceded when input
#'    `GRanges` has 'feature_type' column with proper numeric ranges,
#'    for example 'Promoter (-2000,+200)' or 'TSS(-1000,+1000)'.
#' @param positive `character` vector indicating `GenomicRanges::strand()`
#'    values processed as positive strands, default assumes
#'    either '+' or '*' are positive strand. There should
#'    be no features with '*'.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' gr <- GenomicRanges::GRanges(
#'    seqnames=c("chr1", "chr1", "chr1", "chr2", "chr2",
#'       "chr3", "chr3", "chr3", "chr4"),
#'    ranges=IRanges::IRanges(
#'       start=c(101, 1, 1, 1,
#'          51, 1, 10001, 10001, 101),
#'       end=c(2100, 1600, 1600, 2200,
#'          480, 2000, 11600, 11600, 2100)),
#'    strand=c("+", "+", "-", "+", "+", "-", "+", "-", "-"),
#'    feature_type=c(
#'       "Promoters (-1500,+500)",
#'       "Promoters (-1500,+500)",
#'       "Promoters (-1500,+500)",
#'       "Promoters (-2000,+200)",
#'       "exons",
#'       "Promoters (-1500,+500)",
#'       "Promoters (-1500,+500)",
#'       "Promoters (-1500,+500)",
#'       "TTS (-1000,+1000)"
#'    )
#' )
#' gr2 <- promoter_to_tss(gr)
#' gr2
#' 
#' @export
promoter_to_tss <- function(
   gr,
   use_feature_type=TRUE,
   upstream=2000,
   downstream=200,
   minimum_width=NULL,
   tts_upstream=1000,
   tts_downstream=1000,
   positive=c("+", "*"),
   verbose=FALSE,
   ...
) {
   # assume all peaks are "Promoter" until we see otherwise
   use_type <- rep("Promoter", length(gr));
   k <- seq_along(gr);
   use_upstream <- rep_len(upstream, length(gr))
   use_downstream <- rep_len(downstream, length(gr))

   # Check for 'feature_type' column
   if (isTRUE(use_feature_type) &&
      "feature_type" %in% colnames(GenomicRanges::values(gr))) {
      use_upstream <- rep_len(0, length(gr))
      use_downstream <- rep_len(0, length(gr))
      use_type <- ifelse(
         grepl(
            "TTS",
            ignore.case=TRUE,
            GenomicRanges::values(gr)$feature_type),
         "TTS",
         ifelse(
            grepl(
               "Promoter",
               ignore.case=TRUE,
               GenomicRanges::values(gr)$feature_type),
            "Promoter",
            "Other"
         )
      )
      k <- jamba::igrep("Promoter|TTS", use_type)
      # use distances
      uft <- unique(GenomicRanges::values(gr[k])$feature_type)
      use_dist <- gsub("^.+[(]|[)][ ]*$|[-+]", "", uft)
      names(use_dist) <- uft;
      if (verbose) {
         jamba::printDebug("promoter_to_tss(): ",
            "use_dist:", use_dist);
      }
      use_dist_m <- data.frame(jamba::rbindList(
         strsplit(unique(use_dist), ","),
         newColnames=c("up", "down")
      ))
      rownames(use_dist_m) <- names(use_dist)
      use_dist_m[,1] <- as.numeric(use_dist_m[,1])
      use_dist_m[,2] <- as.numeric(use_dist_m[,2])
      if (verbose) {
         jamba::printDebug("promoter_to_tss(): ",
            "use_dist_m:");
         print(use_dist_m)
      }
      for (irow in rownames(use_dist_m)) {
         k1 <- (GenomicRanges::values(gr)$feature_type %in% irow);
         use_upstream[k1] <- use_dist_m[irow, 1]
         use_downstream[k1] <- use_dist_m[irow, 2]
      }
   }
   if (verbose) {
      jamba::printDebug("promoter_to_tss(): ",
         "table(use_upstream):");
      print(table(use_upstream));
      jamba::printDebug("promoter_to_tss(): ",
         "table(use_downstream):");
      print(table(use_downstream));
   }
   ## Purpose is to convert coordinates of promoters back to original TSS
   if (length(minimum_width) == 0) {
      minimum_width <- (use_upstream + use_downstream)
   }
   # process non-zero width peaks large enough for the definition
   to_process_all <- (minimum_width > 0);
   to_process <- (
      minimum_width > 0 &
      GenomicRanges::width(gr) >= minimum_width);
   non_process <- (to_process_all & !to_process)
   if (verbose) {
      jamba::printDebug("promoter_to_tss(): ",
         "to_process_all:", to_process_all);
      jamba::printDebug("promoter_to_tss(): ",
         "to_process:", to_process);
      jamba::printDebug("promoter_to_tss(): ",
         "non_process:", non_process);
      # jamba::printDebug("promoter_to_tss(): ",
      #    "table(to_process_all, to_process):");
      # print(table(to_process_all=to_process_all, to_process=to_process));
   }

   # Process entries with at least the minimum width, on "+" or "-" strand
   gr2pos <- (as.character(GenomicRanges::strand(gr)) %in% positive);
   newstart <- ifelse(gr2pos[to_process],
      GenomicRanges::start(gr[to_process]) + use_upstream[to_process],
      GenomicRanges::start(gr[to_process]) + use_downstream[to_process] - 1);
   newend <- ifelse(gr2pos[to_process],
      GenomicRanges::end(gr[to_process]) - use_downstream[to_process] + 1,
      GenomicRanges::end(gr[to_process]) - use_upstream[to_process]);
   if (verbose) {
      jamba::printDebug("promoter_to_tss(): ",
         "gr2pos:", gr2pos);
      print(data.frame(newstart, newend));
   }
   GenomicRanges::ranges(gr[to_process]) <- IRanges::IRanges(
      start=newstart,
      end=newend);
   
   # Handle entries too small for processing above
   if (any(non_process)) {
      # front_clip: The region start is at 0 or 1 of the chromosome,
      # which means the size was truncated at the chromosome boundary.
      front_clip <- GenomicRanges::start(gr[non_process]) <= 1;
      # A front_clip happens when one end of the feature is at
      # the front of a chromosome, suggesting it was clipped here.
      # Otherwise assume it was clipped at the far end of the chromosome.
      # For front_clip trust the end(),
      # otherwise trust the start().
      if (verbose) {
         jamba::printDebug("promoter_to_tss(): ",
            "front_clip:", front_clip);
         print(data.frame(
            start=GenomicRanges::start(gr[non_process]),
            end=GenomicRanges::end(gr[non_process]),
            strand=GenomicRanges::strand(gr[non_process]),
            gr2pos=gr2pos[non_process],
            front_clip=front_clip,
            use_upstream[non_process],
            use_downstream[non_process]
         ))
      }
      newstart_small <- ifelse(front_clip,
         ifelse(gr2pos[non_process],
            GenomicRanges::end(gr[non_process]) - use_downstream[non_process] + 1,
            GenomicRanges::end(gr[non_process]) - use_upstream[non_process] + 1),
         ifelse(gr2pos[non_process],
            GenomicRanges::start(gr[non_process]) + use_upstream[non_process],
            GenomicRanges::start(gr[non_process]) + use_downstream[non_process] + 0
         )
      )
      newend_small <- ifelse(front_clip,
         ifelse(gr2pos[non_process],
            (GenomicRanges::end(gr[non_process]) - use_downstream[non_process]) + 1,
            (GenomicRanges::end(gr[non_process]) - use_upstream[non_process]) + 1),
         ifelse(gr2pos[non_process],
            GenomicRanges::start(gr[non_process]) + use_upstream[non_process] + 0,
            GenomicRanges::start(gr[non_process]) + use_downstream[non_process] + 0
         )
      );
      if (verbose) {
         jamba::printDebug("promoter_to_tss(): ",
            "new_df:");
         print(data.frame(newstart_small, newend_small));
      }
      GenomicRanges::ranges(gr[non_process]) <- IRanges::IRanges(
         start=newstart_small,
         end=newend_small);
   }
   return(gr);
}
