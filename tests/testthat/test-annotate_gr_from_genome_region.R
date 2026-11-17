test_that("annotate_gr_by_genome_region", {
   # elaborate test data
   genome_regions1 <- GenomicRanges::GRanges(
      seqnames=c("chr1"),
      ranges=IRanges::IRanges(
         start=c(1501, 3001, 3901, 5501, 6001, 6501, 6001),
         end = c(3500, 3900, 5500, 6000, 6500, 7000, 8000),
      ),
      strand=c("+"),
      gene_name="TESTGENE1",
      gene_id="ENSG0000001.1",
      feature_type=c("Promoters (-1500,+500)", "exon", "intron",
         "exon", "intron", "exon", "TTS (-1000,+1000)"
      )
   )
   genome_regions2 <- genome_regions1;
   GenomicRanges::end(genome_regions2) <- GenomicRanges::end(genome_regions1) + 35000;
   GenomicRanges::start(genome_regions2) <- GenomicRanges::start(genome_regions1) + 35000;
   GenomicRanges::values(genome_regions2)$gene_name <- "TESTGENE2";
   GenomicRanges::values(genome_regions2)$gene_id <- "ENSG0000002.1"

   start3 <- 28001 - GenomicRanges::end(genome_regions1);
   end3 <- 28001 - GenomicRanges::start(genome_regions1);
   genome_regions3 <- GenomicRanges::GRanges(
      seqnames=c("chr1"),
      ranges=IRanges::IRanges(start=rev(start3), end=rev(end3)),
      strand="-",
      gene_name="TESTGENE3",
      gene_id="ENSG0000003.1",
      feature_type=rev(
         GenomicRanges::values(genome_regions1)$feature_type)
   )

   genome_regions <- c(genome_regions1, genome_regions2, genome_regions3)

   peaks <- GenomicRanges::GRanges(
      seqnames=c("chr1"),
      ranges=IRanges::IRanges(
         start=c(1201, 1801, 3841, 4501, 10001, 19001, 30001, 36001),
         end = c(1300, 1900, 3940, 4600, 10100, 19100, 30100, 36100),
      ),
      strand=c("+"),
      gene_name=paste0("peak_00", 1:8),
      gene_id=paste0("peaks"),
      feature_type="peak"
   )
   peaks_anno <- annotate_gr_by_genome_region(peaks, genome_regions)
   peaks_anno
   testthat::expect_equal(
      peaks_anno$nearest_tss_gene_distance,
      c(-1700, -1100, 839, 1499,
      6999, 5899, -5000, -1900)
   )
   testthat::expect_equal(
      peaks_anno$nearest_gene_distance,
      c(-1700, -1100, 0, 0,
      3000, 1900, -5000, -1900)
   )
   testthat::expect_equal(
      peaks_anno$nearest_gene_name,
      c("TESTGENE1", "TESTGENE2", "TESTGENE3")[c(1, 1, 1, 1, 1, 3, 3, 2)]
   )
   testthat::expect_equal(
      peaks_anno$gene_name_winner,
      c(NA, "TESTGENE1")[c(1, 2, 2, 2, 1, 1, 1, 1)]
   )
   testthat::expect_equal(
      peaks_anno$feature_type_winner,
      c("extragenic", "Promoters (-1500,+500)", "exon",
      "intron")[c(1, 2, 3, 4, 1, 1, 1, 1)]
   )
   testthat::expect_equal(
      peaks_anno$gene_feature_type,
      c(NA, "TESTGENE1 Promoters (-1500,+500)",
      "TESTGENE1 exon,TESTGENE1 intron",
      "TESTGENE1 intron")[c(1, 2, 3, 4, 1, 1, 1, 1)]
   )
   testthat::expect_equal(
      peaks_anno$feature_type,
      c("extragenic",
      "Promoters (-1500,+500)",
      "exon,intron",
      "intron")[c(1, 2, 3, 4, 1, 1, 1, 1)]
   )

})
