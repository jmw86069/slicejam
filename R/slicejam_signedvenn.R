
#' Signed Venn overlap logic
#' 
#' Signed Venn overlap logic
#' 
#' This function is deprecated in favor of `venndir::venndir()`
#' and `venndir::signed_overlaps()`.
#' 
#' @family slicejam deprecated
#' 
#' @export
signed_venn <- function
(setlist,
   overlap_type=c("any", "each", "concordance", "concordant"),
   return_type=c("counts", "list"),
   ...)
{
   ##
   overlap_type <- match.arg(overlap_type);
   
   ## 1sec
   svims <- list2im_signed(setlist, do_sparse=TRUE);

   ## 0.02sec
   if (is.character(svims[1,1])) {
      svimsl <- (svims != "") * 1
      svimsl[is.na(svims)] <- 0;
      svimss <- do.call(paste, lapply(seq_len(ncol(svims)), function(i){
         jamba::rmNA(naValue="", svims[,i]);
      }))
   } else {
      svimsl <- (svims != 0) * 1
      svimsl[is.na(svims)] <- 0;
      svimss <- do.call(paste, lapply(seq_len(ncol(svims)), function(i){
         jamba::rmNA(naValue="0", svims[,i]);
      }))
   }

   ## 1.2sec
   svimsv <- do.call(paste, lapply(seq_len(ncol(svimsl)), function(i){svimsl[,i]}))
   
   ## concordance
   ## 0.02sec
   svimssu <- unique(svimss);
   svimssu_concordance <- jamba::nameVector(sapply(strsplit(svimssu, " "), function(i){
      j <- setdiff(i, c("","0"))
      length(unique(j)) == 1;
   }), svimssu)

   ## split by observed directions within each overlap set
   ## 0.07 sec
   svims_split <- split(svimss, svimsv);
   
   ## Create labels
   svims_split_names <- sapply(jamba::nameVectorN(svims_split), function(i){
      paste(collapse="&",
         names(setlist)[strsplit(i, " ")[[1]] %in% "1"])
   })
   svims_df <- data.frame(jamba::rbindList(lapply(strsplit(names(svims_split), " "), as.numeric)),
      check.names=FALSE);
   colnames(svims_df) <- names(setlist);
   svims_df$sum <- rowSums(svims_df[,names(setlist),drop=FALSE]);
   svims_df <- jamba::mixedSortDF(svims_df, byCols=c("sum", paste0("-", names(setlist))))
   rownames(svims_df) <- apply(svims_df[,names(setlist),drop=FALSE], 1, function(i){
      paste(collapse="&",
         names(setlist)[i != 0])
   })
   svims_split_names1 <- jamba::nameVector(rownames(svims_df),
      jamba::pasteByRow(svims_df[,names(setlist),drop=FALSE], sep=" ", condenseBlanks=FALSE));
   svims_split_names <- jamba::nameVector(
      jamba::pasteByRow(svims_df[,names(setlist),drop=FALSE], sep=" ", condenseBlanks=FALSE),
      rownames(svims_df))
   
   ret_vals <- list();
   
   # 0.08sec
   if ("counts" %in% return_type) {
      svims_split_counts <- lapply(svims_split_names, function(iname){
         i <- svims_split[[iname]];
         if ("concordance" %in% overlap_type) {
            j <- ifelse(svimssu_concordance[i], i, "mixed")
         } else if ("concordant" %in% overlap_type) {
            j <- ifelse(svimssu_concordance[i], "concordant", "mixed")
         } else if ("any" %in% overlap_type) {
            j <- rep(iname, length(i));
         } else {
            j <- i;
         }
         #jv <- table(j);
         #nameVector(as.vector(jv), names(jv))
         #data.frame(count=as(rev(table(j)), "matrix"), check.names=FALSE)
         jdf <- data.frame(count=as(table(j), "matrix"), check.names=FALSE);
         if (nrow(jdf) > 1) {
            jrows <- rev(jamba::provigrep(c("mixed", "concordan", "."),
               jamba::mixedSort(rownames(jdf))));
            jdf[jrows,,drop=FALSE];
         } else {
            jdf;
         }
      });
      ret_vals$Venn_Counts <- svims_split_counts;
   }
   
   ## list of list of vectors
   if ("list" %in% return_type) {
      svims_split_list <- lapply(svims_split_names, function(iname){
         i <- svims_split[[iname]];
         if ("concordance" %in% overlap_type) {
            j <- ifelse(svimssu_concordance[i], i, "mixed")
         } else if ("concordant" %in% overlap_type) {
            j <- ifelse(svimssu_concordance[i], "concordant", "mixed")
         } else if ("any" %in% overlap_type) {
            j <- rep(iname, length(i));
         } else {
            j <- i;
         }
         #rev(split(i, j));
         split(i, j);
      });
      ret_vals$Venn_List <- svims_split_list;
   }
   return(ret_vals);
}

#' Label Venn polygons
#' 
#' Label Venn polygons
#' 
#' This function is deprecated in favor of equivalent functions
#' in the R package `venndir`, such as `venndir::polygon_label_outside()`.
#' 
#' @family slicejam deprecated
#' 
#' @export
label_venn_polys <- function
(venn_spdf,
 venn_counts=NULL,
 label_style=1,
 ...)
{
   vps <- gridBase::baseViewports();
   grid::pushViewport(vps$inner, vps$figure, vps$plot);
   on.exit(grid::popViewport(3));

   ## Few different label style options
   #
   # 1 - colored labels with border
   up1 <- "#990000";
   dn1 <- "#000099";
   up2 <- "#ff8888";
   dn2 <- "#8888ff";
   label_n <- length(venn_spdf$color);
   if (1 %in% label_style) {
      label_fill <- venn_spdf$color;
      label_border <- jamba::makeColorDarker(venn_spdf$color,
         darkFactor=1.2);
      label_color <- jamba::setTextContrastColor(venn_spdf$color);
   }
   # - colored shaded labels without border
   if (2 %in% label_style) {
      label_fill <- jamba::alpha2col(alpha=0.8,
         venn_spdf$color);
      label_border <- NA;
      label_color <- jamba::setTextContrastColor(venn_spdf$color);
   }
   # - light background labels with border
   if (3 %in% label_style) {
      label_fill <- rep("#ffeeaa", length.out=label_n);
      label_border <- jamba::alpha2col(alpha=0.8,
         jamba::makeColorDarker(venn_spdf$color,
            darkFactor=1.2));
      label_color <- rep("#000000", length.out=label_n);
      up2 <- up1;
      dn2 <- dn1;
   }
   # - light shaded background labels without border
   if (4 %in% label_style) {
      label_fill <- rep("#ffeeaa99", length.out=label_n);
      label_border <- NA;
      label_color <- rep("#000000", length.out=label_n);
      up2 <- up1;
      dn2 <- dn1;
   }

   ## glyph lookup table
   ## - convert "concordant" to a symbol to display
   ## - convert "mixed"
   
   
   ## venn_counts
   if (length(venn_counts) == 0) {
      venn_counts_base <- venn_spdf$venn_counts;
   } else {
      venn_counts_base <- sapply(venn_counts[venn_spdf1$label], function(idf){
         sum(idf$count)
      });
   }
   
   ## Central text label
   venn_text <- ifelse(grepl("&", venn_spdf$label),
      jamba::formatInt(jamba::rmNA(naValue=0,
         venn_counts_base)),
      paste0("**",
         venn_spdf$label,
         "**<br>\n",
         jamba::formatInt(jamba::rmNA(naValue=0,
            venn_counts_base))));
   g_text <- gridtext::richtext_grob(
      text=venn_text,
      x=venn_spdf$x_label,
      y=venn_spdf$y_label,
      default.units="native",
      vjust=0.5,
      hjust=1,#ifelse(nchar(up_text) == 0 & nchar(dn_text) == 0, 0.5, 1),
      halign=0.5,
      rot=0,
      padding=grid::unit(c(2), "pt"),
      r=grid::unit(c(2), "pt"),
      gp=grid::gpar(col=label_color),
      box_gp=grid::gpar(
         col=label_border,
         fill=label_fill,
         lty=1)
   );
   grid::grid.draw(g_text)

   ## prepare up and down arrow labels
   # symbol2utf8("upArrow")
   up1 <- "#990000";
   dn1 <- "#000099";
   up2 <- "#ff8888";
   dn2 <- "#8888ff";
   upArrow <- "\u2191";
   dnArrow <- "\u2193";
   doubleUpArrow <- "\u21c8";
   doubleDnArrow <- "\u21ca";
   upDnArrow <- "\u21c5";

   parallelTo <- "\u2191\u2191";
   notParallelTo <- "\u2191\u2193";
   #parallelTo <- "\u21c9";
   #notParallelTo <- "\u21b9";
   
   upArrow_text_b <- paste0("<span style='color:", up1, "'>**",
      upArrow,
      "**</span>");
   upArrow_text_d <- paste0("<span style='color:", up2, "'>**",
      upArrow,
      "**</span>");
   conArrow_text_b <- paste0("<span style='color:", dn1, "'>**",
      parallelTo,
      "**</span>");
   conArrow_text_d <- paste0("<span style='color:", dn2, "'>**",
      parallelTo,
      "**</span>");
   
   dnArrow_text_b <- paste0("<span style='color:", dn1, "'>**",
      dnArrow,
      "**</span>");
   dnArrow_text_d <- paste0("<span style='color:", dn2, "'>**",
      dnArrow,
      "**</span>");
   disArrow_text_b <- paste0("<span style='color:", up1, "'>**",
      notParallelTo,
      "**</span>");
   disArrow_text_d <- paste0("<span style='color:", up2, "'>**",
      notParallelTo,
      "**</span>");
   venn_bright <- (jamba::col2hcl(venn_spdf$color)["L",] >= 55);
   
   ## Check for any entries with more than one category
   if (any(sdim(venn_counts)$rows > 1)) {
      extra_texts <- sapply(jamba::nameVectorN(venn_counts), function(iname){
         ivenn <- venn_counts[[iname]];
         if (!grepl("[&]", iname) & grepl("concordan", rownames(ivenn))) {
            ## Skip "concordant" label for single entry
            NULL;
         } else {
            glyph_text <- ifelse(grepl("concordan", rownames(ivenn)),
               conArrow_text_b,
               ifelse(grepl("mixed", rownames(ivenn)),
                  disArrow_text_b,
                  ""));
            label_texts <- ifelse(
               jamba::rmNA(naValue=0, ivenn$count) > 0,
               paste(
                  ifelse(venn_bright,
                     glyph_text,
                     glyph_text),
                  jamba::rmNA(naValue=0, formatInt(ivenn$count))),
               "");
            paste(label_texts, collapse="<br>");
         }
      })
      has_extra <- which(sapply(extra_texts, nchar) > 0);
      if (any(has_extra)) {
         extra_texts_use <- unlist(extra_texts[has_extra]);
         g_extra <- gridtext::richtext_grob(
            text=extra_texts_use,
            x=unlist(venn_spdf$x_label[has_extra]),
            y=unlist(venn_spdf$y_label[has_extra]),
            default.units="native",
            vjust=0.5,
            hjust=0,
            halign=0,
            rot=0,
            padding=grid::unit(c(2), "pt"),
            r=grid::unit(c(2), "pt"),
            gp=grid::gpar(
               col=label_color[has_extra],
               fontsize=8,
               fontfamily="Arial"),
            box_gp=grid::gpar(
               col=label_border[has_extra],
               fill=label_fill[has_extra],
               lty=1)
         )
         grid::grid.draw(g_extra);
      }
   }

   ## Close the grid viewport on this base plot panel
   #popViewport(3);
   
   return(invisible(NULL));
   
   up_text <- ifelse(
      floor(jamba::rmNA(naValue=0, venn_spdf$venn_counts)/20) > 0,
      paste(
         ifelse(venn_bright,
            upArrow_text_b,
            upArrow_text_d),
         floor(venn_spdf$venn_counts/20)),
      "");
   dn_text <- ifelse(
      floor(jamba::rmNA(naValue=0, venn_spdf$venn_counts)/25) > 0,
      paste(
         ifelse(venn_bright,
            dnArrow_text_b,
            dnArrow_text_d),
         floor(venn_spdf$venn_counts/25)),
      "");

   ## Up labels
   has_up <- nchar(up_text) > 0;
   g_up <- gridtext::richtext_grob(
      text=up_text[has_up],
      x=venn_spdf$x_label[has_up],
      y=venn_spdf$y_label[has_up],
      default.units="native",
      vjust=0,
      hjust=0,
      halign=0,
      rot=0,
      padding=grid::unit(c(2), "pt"),
      r=grid::unit(c(2), "pt"),
      gp=grid::gpar(
         col=label_color[has_up],
         fontsize=12),
      box_gp=grid::gpar(
         col=label_border[has_up],
         fill=label_fill[has_up],
         lty=1)
   )
   
   ## Down labels
   has_dn <- nchar(dn_text) > 0;
   g_dn <- gridtext::richtext_grob(
      text=dn_text[has_dn],
      x=venn_spdf$x_label[has_dn],
      y=venn_spdf$y_label[has_dn],
      default.units="native",
      vjust=1,
      hjust=0,
      halign=0,
      rot=0,
      padding=grid::unit(c(2), "pt"),
      r=grid::unit(c(2), "pt"),
      gp=grid::gpar(
         col=label_color[has_dn],
         fontsize=12),
      box_gp=grid::gpar(
         col=label_border[has_dn],
         fill=label_fill[has_dn],
         lty=1)
   )
   
   ## Close the grid viewport on this base plot panel
   popViewport(3);
   
}

