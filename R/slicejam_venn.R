
## slicejam_venn.R

#' Directional Venn diagram
#' 
#' Directional Venn diagram
#' 
#' @export
venndir <- function
(iHitList,
   fill_colors=NULL,
   do_signed=TRUE,
   do_labels=FALSE,
   do_proportional=FALSE,
   do_arrows_only=FALSE,
   non_overlapping_labels=FALSE,
   srt_labels=-20,
   label_poly_scale=0.92,
   label_cex=NULL,
   main="Venn Diagram",
   mtext_line=-4,
   mtext_quantile_pos=0.4,
   sub_color=c("grey20", "grey20"),
   sub_font=c(2,1),
   dn_order=c(2,1,3),
   remove_grep="^A_[0-9]+_|^LOC[0-9]{4,}|^DKFZ[0-9]{2,}|^FLJ[0-9]{4,}|^kiaa[0-9]|^c.*orf",
   use_subset=NULL,
   subset_name="",
   return_top_entries=10,
   do_caption=TRUE,
   verbose=FALSE,
   ...)
{
   ## Purpose is a simple wrapper to draw signed Venn diagrams.
   ## If statsHitsA is a list object, with an element named "statsHitsA" then
   ## that element will be used for the Venn.  This option allows using the
   ## output of se_contrast_stats() directly.
   ##
   ## doSigned=TRUE will use signedVennPlotSets()
   ## doSigned=FALSE will only use the names and then call vennPlotSets()
   ##
   ## doLabels=TRUE will call vennPlotSets2(vennPlotSets()) on the names
   ##
   ## doProportional=TRUE will not include labels within the Venn diagram,
   ## but will call plotVenneuler() which creates proportional Venn diagrams
   ## (Euler diagrams). Note that they do not work very well with 4 or more sets.#
   ##
   ## returnTopEntries>0 will return the top 'n' entries per Venn subset, but only
   ## if supplying allNormStats and not allNormStats$statsHitsA.  It will use
   ## the P-value and fold change for ordering the list as relevant.
   
   mtextQuantilePos <- rep(mtextQuantilePos, length.out=2);
   
   allNormStats <- NULL;
   
   ## Optional shortcut: accept simply the iHitList, similar to sending
   ## a signed setlist to signedVennPlotSets()
   

   ## Optionally remove certain gene patterns
   if (length(remove_grep) > 0) {
      if (verbose) {
         jamba::printDebug("venndir(): ",
            "remove_grep:",
            paste0("'", remove_grep, "'"));
      }
      iHitList <- lapply(jamba::nameVectorN(iHitList), function(i_name){
         i <- iHitList[[i_name]];
         i_names <- names(i);
         drop_names <- jamba::provigrep(remove_grep,
            i_names);
         i_keep <- i[!names(i) %in% drop_names];
         if (verbose && length(drop_names) > 0) {
            jamba::printDebug("venndir(): ",
               "remove_grep removed ",
               jamba::formatInt(length(drop_names)),
               " of ",
               jamba::formatInt(length(i)),
               " from ",
               i_name);
         }
         i_keep;
      });
   }
   
   ########################################################################
   ## Optionally only include a subset of genes provided by iGenesSubset
   if (length(use_subset) > 0) {
      if (verbose) {
         jamba::printDebug("venndir(): ",
            "Applying use_subset");
      }
      iHitList <- lapply(jamba::nameVectorN(iHitList), function(i_name){
         i <- iHitList[[i_name]];
         i_names <- names(i);
         i_drop <- (!i_names %in% use_subset);
         iKeep <- intersect(iNamesG1, iGenesSubset);
         if (verbose && any(i_drop)) {
            jamba::printDebug("venndir(): ",
               "use_subset removed ",
               jamba::formatInt(sum(i_drop)),
               " of ",
               jamba::formatInt(length(i)),
               " from ",
               i_name);
         }
         i[iKeep];
      });
      iSubmain1nl <- length(strsplit(iSubmain1, "\n")[[1]]);
      iSubmain2nl <- length(strsplit(iSubmain2, "\n")[[1]]);
      iSubmain1nlsep <- "\n";
      iSubmain2nlsep <- "\n";
      if (iSubmain1nl > iSubmain2nl) {
         iSubmain1nlsep <- paste0(iSubmain2nlsep,
            rep("\n", length.out=(iSubmain1nl-iSubmain2nl)), collapse="");
      } else if (iSubmain1nl < iSubmain2nl) {
         iSubmain2nlsep <- paste0(iSubmain1nlsep,
            rep("\n", length.out=(iSubmain2nl-iSubmain1nl)), collapse="");
      }
      iSubmain1 <- paste0(iSubmain1, iSubmain1nlsep, "Subset:");
      subsetLabel <- paste0(iGenesSubsetName,
         "(", formatInt(length(iGenesSubset)), " genes)");
      iSubmain2 <- paste0(iSubmain2, iSubmain2nlsep, subsetLabel);
   }
   
   ## Return data in form of a list
   retList <- list();
   retList$inputHitList <- iHitList;
   
   ## Define colors
   if (length(fill_colors) < length(iHitList)) {
      fill_colors <- colorjam::rainbowJam(length(iHitList));
   }
   
   ## Proportional Venn (Euler) diagrams currently do not support
   ## using labels
   if (doProportional) {
      if (!suppressPackageStartupMessages(require(venneuler))) {
         do_proportional <- FALSE;
         jamba::printDebug("venndir(): ",
            c("Warning: The '",
               "venneuler",
            "' package is required for ",
            "doProportional=TRUE.",
            " Proceeding with ", 
               "doProportional=FALSE."),
            sep="");
      }
      do_labels <- FALSE;
      do_signed <- FALSE;
   }
   
   if (!do_signed || do_labels) {
      ## Special case with signed Venn and labels, we will create a prefix
      ## with up/down arrows to each entry
      if (do_signed && do_labels) {
         ## Signed, use the names, but only after putting arrows into the labels
         iHitList <- signed_venn_counts(iHitList,
            returnType="setlist_arrows",
            ...);
         #names(iHitList) <- gsub("!", "-", names(iHitList));
         #retList$iHitList <- iHitList;
      } else {
         ## Not directional, use only the names
         iHitList <- lapply(iHitList, names);
      }
      if (do_labels) {
         ## Make Venn diagram, then make a Venn diagram with labels inside the shapes
         if (verbose) {
            jamba::printDebug("venndir(): ",
               "Make a Venn diagram with labels inside the shapes.");
         }
         vps1 <- venn_plot_sets(iHitList,
            mymain="",
            verbose=verbose,
            doPlot=FALSE);
         if (verbose) {
            printDebug("venndir(): ",
               "Running venn_plot_sets2 with labels.");
         }
         vps <- venn_plot_sets2(vps1,
            labelColor=attr(iHitList, "color"),
            verbose=verbose,
            doNonOverlappingLabels=doNonOverlappingLabels,
            srt=srtLabels,
            labelCex=labelCex,
            fill_colors=fill_colors,
            labelPolyScale=labelPolyScale,
            doArrowsOnly=doArrowsOnly,
            ...);
      } else if (do_proportional) {
         ## Do Venn diagram, then a proportional Venn diagram (Euler diagram)
         if (verbose) {
            jamba::printDebug("venndir(): ",
               "Do Venn diagram, then a proportional Venn diagram (Euler diagram)");
         }
         vps <- venn_plot_sets(iHitList,
            mymain="",
            doPlot=FALSE,
            ...);
         vps1 <- plot_venneuler(iHitList,
            mymain="",
            ...);
         vps$vennCircles <- vps1$ellList;
      } else {
         ## Do normal Venn diagram, do not put labels into the plot
         if (verbose) {
            jamba::printDebug("venndir(): ",
               "Do normal Venn diagram, do not put labels into the plot");
         }
         vps <- venn_plot_sets(iHitList,
            mymain="",
            lcol=fillColors,
            lines=fillColors,
            ...);
      }
   } else {
      ## Do signed venn diagram, no labels inside
      if (verbose) {
         jamba::printDebug("venndir(): ",
            "Do signed venn diagram, no labels inside.");
      }
      vps <- signed_venn_plot_sets(iHitList,
         mymain="",
         lcol=fillColors,
         lines=fillColors,
         ...);
      imSigned <- list2im_signed(iHitList);
      imFactor <- jamba::pasteByRow(imSigned, sep=":");
      if (verbose) {
         jamba::printDebug("venndir(): ",
            "Completed signed venn diagram, no labels inside.");
      }
   }
   title(main=iMain,
      font.main=1,
      line=2,
      cex.main=1.75);
   if (do_caption) {
      mtext(line=mtextLine,
         at=quantile(par("usr")[1:2], mtextQuantilePos),
         c(iSubmain1,
            iSubmain2),
         col=subMainColor,
         adj=c(1,0),
         font=subMainFont);
   }
   retList$vps <- vps;
   retList$iHitList <- iHitList;
   invisible(retList);
}

#' Define directional Venn counts
#' 
#' Define directional Venn counts
#' 
#' @examples
#' L1 <- list(A=c("C","A","B","A"),
#'    D=c("D","E","F","D"),
#'    A123=LETTERS[c(1:8,3,5)],
#'    T=LETTERS[7:9]);
#' L1;
#' # Convert each vector to a signed vector
#' set.seed(123);
#' L2 <- lapply(L1, function(i){
#'    i <- unique(i);
#'    jamba::nameVector(sample(c(-1,1), size=length(i), replace=TRUE), i);
#' });
#' L2;
#' 
#' signed_venn_counts(L2)
#' signed_venn_counts(L2, return_type="setlist_arrows")
#' 
#' @export
signed_venn_counts <- function
(setlist,
 return_type=c("counts", "names", "im", "setlist_arrows",
    "dircounts"),
 sep="-",
 alt_sep="!",
 pos_col="red4",
 neg_col="blue4",
 mixed_col=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take a list of named numeric vectors and return the overlap counts,
   ## including the counts of all-up, all-down, and mixed directionality
   ##
   ## returnType="setlistWithArrows" will return setlist except prepending UTF8 arrows so
   ## the elements will indicate directionality of the original data. This setlist can be
   ## used in standard Venn diagrams which label the elements.
   return_type <- match.arg(return_type);
   ## Make sure names are unique
   names(setlist) <- jamba::makeNames(names(setlist));
   
   ## Correct names which may include the separator to be used later on
   names(setlist) <- gsub(sep,
      alt_sep,
      names(setlist),
      fixed=TRUE);
   if (verbose) {
      jamba::printDebug("signed_venn_counts(): ",
         "list2im_signed()");
   }
   setlistim <- list2im_signed(setlist);
   if (return_type %in% "im") {
      return(setlistim);
   }
   if (return_type %in% "setlist_arrows") {
      if (verbose) {
         jamba::printDebug("signed_venn_counts(): ",
            "setlist_arrows");
      }
      upArrow <- symbol2utf8("upArrow");
      downArrow <- symbol2utf8("downArrow");
      
      ## Create unique row matrix subset to improve efficiency
      setlistim_u <- unique(setlistim);
      u_match <- match.im(setlistim, setlistim_u);
      if (verbose) {
         jamba::printDebug("signed_venn_counts(): ",
            "Completed u_match, length(u_match):",
            length(u_match),
            ", nrow(setlistim_u):",
            nrow(setlistim_u));
      }
      
      ## Now define colors
      setlistim_neg <- rowSums(setlistim_u < 0);
      setlistim_pos <- rowSums(setlistim_u > 0);

      ## Set of vector of directional colors, then split by name,
      ## then blend colors
      setlistim_pos_colv <- rep(rep(pos_col, length(setlistim_pos)), setlistim_pos);
      setlistim_pos_coln <- rep(seq_along(setlistim_pos), setlistim_pos);
      setlistim_neg_colv <- rep(rep(neg_col, length(setlistim_neg)), setlistim_neg);
      setlistim_neg_coln <- rep(seq_along(setlistim_neg), setlistim_neg);
      setlistim_col_l <- split(
         c(setlistim_pos_colv, setlistim_neg_colv),
         c(setlistim_pos_coln, setlistim_neg_coln));
      ## Blend colors
      element_colors_u <- blend_colors(setlistim_col_l);
      element_colors <- element_colors_u[u_match];
      #setlistim_pos_neg_col_bl_u <- blend_colors(setlistim_col_l);
      #setlistim_pos_neg_col_bl <- setlistim_pos_neg_col_bl_u[u_match];
      if (verbose) {
         jamba::printDebug("signed_venn_counts(): ",
            "Completed color blend per entry, length(setlistim_col_l):",
            length(setlistim_col_l));
         
      }

      ## Optionally override some mixed colors
      if (length(mixed_col) > 0 && any(!setlistim_pos_neg_col_bl %in% c(pos_col, neg_col))) {
         #jamba::printDebug("Assigning mixed_col");
         which_mixed <- (!setlistim_pos_neg_col_bl %in% c(pos_col, neg_col));
         setlistim_pos_neg_col_bl[which_mixed] <- mixed_col;
      }
      
      ## Generate vector of up/down arrows for each entry
      ## operate on unique matrix to reduce processing
      element_arrows_u <- apply(setlistim_u, 1, function(i){
         #paste(ifelse(i == -1, downArrow, ifelse(i == 1, upArrow, "-")), collapse="")
         ## Changed 08jul2016 to use sign instead of integer values
         paste(ifelse(i < 0, downArrow, ifelse(i > 0, upArrow, "-")), collapse="");
      });
      element_arrows <- element_arrows_u[u_match];
      if (verbose) {
         jamba::printDebug("signed_venn_counts(): ",
            "Completed element_arrows.");
         jamba::printDebug("length(element_arrows):", length(element_arrows));
         jamba::printDebug("nrow(setlistim):", nrow(setlistim));
      }
      
      new_names <- paste(element_arrows, rownames(setlistim));
      if (1 == 2) {
         setlistin_pos_neg_col_bl_df <- data.frame(stringsAsFactors=FALSE,
            check.names=FALSE,
            name=rownames(setlistim),
            new_name=new_names,
            color=setlistim_pos_neg_col_bl);
         attr(setlist_with_arrows, "color") <- setlistin_pos_neg_col_bl_df;
      }

      rownames(setlistim) <- new_names;
      if (verbose) {
         jamba::printDebug("venndir(): ",
            "Beginning im2list():");
      }
      setlist_with_arrows <- multienrichjam::im2list(setlistim);
      #setlist_with_arrows <- apply(setlistim, 2, function(i){
      #   i_which <- which(i != 0);
      #   paste(element_arrows, rownames(setlistim))[i_which];
      #});
      if (verbose) {
         jamba::printDebug("signed_venn_counts(): ",
            "Completed setlist_with_arrows.");
      }
      ##
      attr(setlist_with_arrows, "element_colors") <- element_colors;
      attr(setlist_with_arrows, "element_arrows") <- element_arrows;
      return(setlist_with_arrows);
   }
   if (verbose) {
      jamba::printDebug("signed_venn_counts(): ",
         "venn_combn");
   }
   
   venn_combn <- unlist(recursive=FALSE, sapply((1:ncol(setlistim)), function(i){
      unlist(recursive=FALSE, apply(combn(colnames(setlistim), i), 2, function(j){
         k <- jamba::nameVector(rep(0, ncol(setlistim)), colnames(setlistim));
         k[j] <- 1;
         ## Note: commented out the rename step below, we do renaming at the top of this function
         #jName <- paste(gsub(sep, altSep, j), collapse="-");
         j_name <- paste(j, collapse=sep);
         k_list <- list(c(k));
         names(k_list) <- j_name;
         k_list;
      }));
   }));
   
   if (verbose) {
      jamba::printDebug("signed_venn_counts(): ",
         "venn_combn_rows");
   }
   i_rowsums <- rowSums(abs(setlistim));
   i_allup <- rowSums(setlistim > 0) == i_rowsums;
   i_alldn <- rowSums(setlistim < 0) == i_rowsums;

   venn_combn_rows <- lapply(venn_combn, function(i_combn){
      i_hit_col <- names(which(i_combn != 0));
      i_non_hit_col <- names(which(i_combn == 0));
      i_hit_rowsums <- rowSums(abs(setlistim[,i_hit_col,drop=FALSE]));
      i_which_hits <- names(which(i_hit_rowsums == length(i_hit_col) &
            i_rowsums == length(i_hit_col)));
      i_which_hits_sub <- setlistim[i_which_hits,i_hit_col,drop=FALSE];
      i_which_hits_all_up <- rownames(i_which_hits_sub)[rowSums(i_which_hits_sub) == ncol(i_which_hits_sub)];
      i_which_hits_all_dn <- rownames(i_which_hits_sub)[rowSums(i_which_hits_sub) == -1*ncol(i_which_hits_sub)];
      i_which_hits_all_mixed <- rownames(i_which_hits_sub)[abs(rowSums(i_which_hits_sub)) != ncol(i_which_hits_sub)];
      list(all_hits=i_which_hits,
         all_up=i_which_hits_all_up,
         all_dn=i_which_hits_all_dn,
         all_mixed=i_which_hits_all_mixed);
   });
   if (return_type %in% "counts") {
      return(lapply(venn_combn_rows, lengths));
   } else {
      return(venn_combn_rows);
   }
}

#' Convert symbol name to UTF8 Unicode
#' 
#' @export
symbol2utf8 <- function
(x,
   ...)
{
   ## Purpose is to wrapper some commonly used symbols to generate UTF8 codes
   ## given only a simple name.
   ##
   ## As a trick, if x=TRUE, it will return the full set of entries defined
   ## by this function.
   ##
   ## NEEDS TO BE EXPANDED
   utf8set <- c(
      "greaterThanEqualTo"="\u2265",
      "lessThanEqualTo"="\u2264",
      "equalTo"="\u003d",
      "notEqualTo"="\u2260",
      "identicalTo"="\u2261",
      "notIdenticalTo"="\u2262",
      "parallelTo"="\u2225",
      "notParallelTo"="\u2226",
      "upArrow"="\u2191",
      "downArrow"="\u2193",
      "leftArrow"="\u2190",
      "rightArrow"="\u2192",
      "doubleUpArrow"="\u21c8",
      "doubleDownArrow"="\u21ca",
      "upDownArrow"="\u21c5",
      "rightLeftArrow"="\u21c4",
      "doubleRightArrow"="\u21c9",
      "male"="\u2642",
      "female"="\u2640");
   if (x %in% c(TRUE)) {
      return(utf8set);
   }
   ## Convert to lowercase, remove spaces, etc.
   x <- nameVector(gsub("[-_. ]+", "", tolower(x)), x);
   if (all(x) %in% tolower(names(utf8set))) {
      y <- utf8set[match(x, tolower(names(utf8set)))];
   } else {
      ## Use grep to match most likely entries
      xGrep <- sapply(x, function(i){
         iGrep <- paste("^", i, sep="");
         iSet <- igrep(iGrep, names(utf8set))[1];
      });
      y <- nameVector(utf8set[xGrep], names(x));
   }
   y;
}

#' Venn plot sets base function
#' 
#' @export
vennPlotSets <- function
(setlist=NULL,
 alpha=0.35,
 showZeros=FALSE,
 mymain="Venn Plot",
 mysub="",
 vennSetList=NULL,
 doPlot=TRUE,
 verbose=FALSE,
 sep="-",
 altSep="!",
 lcol=NULL,
 lines=NULL,
 whichSets=NULL,
 xlim=c(0,10),
 ylim=c(0,10),
 type="ellipse",
 ...)
{
   ## Purpose is to calculate and plot Venn overlaps
   ##
   #if (!suppressPackageStartupMessages(require(sp))) {
   #   stop("vennPlotSets() requires the sp package for Polygon objects.");
   #}
   if (length(setlist) == 0 & length(vennSetList) == 0) {
      return(NULL);
   }
   if (length(setlist) == 0 & length(vennSetList) > 0) {
      # if setlist is missing, assume we're getting vennSetList from a previous run
      setlist <- sapply(colnames(vennSetList$Intersect_Matrix), function(i){
         rownames(vennSetList$Intersect_Matrix)[vennSetList$Intersect_Matrix[,i] %in% c(1)]
      })
   }
   ## Allow alternate input in form of a matrix with 1's and 0's
   if (class(setlist) %in% c("TestResults")) {
      if (verbose) {
         printDebug("vennPlotSets(): ",
            "Converting class ",
            class(setlist),
            " into matrix.");
      }
      setlist <- as(setlist, "matrix");
   }
   if (jamba::igrepHas("matrix|data.frame", class(setlist))) {
      if (verbose) {
         jamba::printDebug("vennPlotSets(): ",
            "Converting class ",
            class(setlist),
            " into list.");
      }
      setlist <- lapply(jamba::nameVector(colnames(setlist)), function(i){
         i1 <- which(!setlist[,i] %in% c(0,"0",FALSE));
         rownames(setlist)[i1];
      });
   }
   if (verbose) {
      jamba::printDebug("vennPlotSets(): ",
         "Beginning overLapper() method for vennsets.");
   }
   
   ## Make sure the list is named, and has unique names
   if (is.null(names(setlist)) || all(is.na(names(setlist)))) {
      names(setlist) <- paste("setlist", 1:length(setlist), sep="_");
   }
   if (max(tcount(names(setlist))) > 1) {
      names(setlist) <- jamba::makeNames(names(setlist));
   }
   ## Define colors here
   ngrp <- length(setlist);
   if (is.null(lines)) {
      lines <- colorjam::rainbowJam(ngrp);
   }
   if (is.null(lcol)) {
      lcol <- jamba::makeColorDarker(darkFactor=1.2, lines);
   }
   
   ## Optionally only evaluate some of the sets
   if (!is.null(whichSets)) {
      setlist <- setlist[whichSets];
      ## accordingly, adjust the colors to be consistent with the
      ## original full list
      lcol <- lcol[whichSets];
      lines <- lines[whichSets];
   }
   
   ## Run the overLapper() method which calculates the set overlaps
   olSets <- overLapper(setlist=setlist,
      type="vennsets",
      verbose=verbose,
      sep=sep,
      ...);
   if (verbose) {
      jamba::printDebug("vennPlotSets(): ",
         "Beginning overLapper() method for intersects.");
   }
   olIntersects <- overLapper(setlist=setlist,
      type="intersects",
      verbose=verbose,
      sep=sep,
      ...);
   if (verbose) {
      jamba::printDebug("vennPlotSets(): ",
         "Generating counts summary.");
   }
   counts <- list(sapply(olSets$Venn_List, length));
   olIntersects$Venn_List <- olSets$Venn_List;
   olIntersects$counts <- counts;
   
   if (verbose) {
      jamba::printDebug("vennPlotSets(): ",
         "Calling vennPlot(..., doPlot=",
         doPlot, ")");
   }
   vp1 <- vennPlot(counts,
      alpha=alpha,
      showZeros=showZeros,
      mymain=mymain,
      mysub=mysub,
      vennList=olSets$Venn_List,
      sep=sep,
      lcol=lcol,
      lines=lines,
      doPlot=doPlot,
      xlim=xlim,
      ylim=ylim,
      type=type,
      ...);
   olIntersects <- c(olIntersects,
      list(vennCircles=vp1,
         xlim=xlim,
         ylim=ylim));

   invisible(olIntersects);
}

#' Compute set overlaps
#' 
#' Compute set overlaps
#' 
#' This function was graciously provided by Dr. Thomas Girke,
#' copied from August 15, 2009 and thereafter adapted by
#' James M. Ward.
#' 
#' @export
overLapper <- function
(setlist=setlist,
 complexity=seq_along(setlist),
 sep="-",
 cleanup=FALSE,
 keepdups=FALSE,
 type=c("vennsets", "intersects"),
 verbose=FALSE,
 ...)
{
   ##########################################
   ## Intersect and Venn Diagram Functions ##
   ##########################################
   ## Author: Thomas Girke
   ## Last update: August 15, 2009
   ## Utilities:
   ## (1) Venn Intersects
   ##     Computation of Venn intersects among 2-20 or more sample sets using the typical
   ##     'only in' intersect logic of Venn comparisons, such as: objects present only in
   ##     set A, objects present only in the intersect of A & B, etc. Due to this restrictive
   ##     intersect logic, the combined Venn sets contain no duplicates.
   ## (2) Regular Intersects
   ##     Computation of regular intersects among 2-20 or more sample sets using the
   ##     following intersect logic: objects present in the intersect of A & B, objects present
   ##     in the intersect of A & B & C, etc. The approach results usually in many duplications
   ##     of objects among the intersect sets.
   ## (3) Graphical Utilities
   ##     - Venn diagrams of 2-5 sample sets.
   ##     - Bar plots for the results of Venn intersect and all intersect approaches derived
   ##       from many samples sets.
   ##
   ## Detailed instructions for using the functions of this script are available on this page:
   ##     http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn
   ##
   ##
   ## Modified by James M. Ward to enable signed set comparisons, using named lists where the class
   ## is either numeric or integer.
   
   #######################################
   ## Define Generic Intersect Function ##
   #######################################
   ## Computation of (1) Venn Intersects and (2) Regular Intersects
   ##
   
   ## Match function arguments
   type <- match.arg(type);
   
   ## First test for signed comparisons
   doSigned <- FALSE;
   if (all(sapply(setlist, class) %in% c("numeric", "integer"))) {
      ## Now check that the lists each are named lists
      namedCheck <- sapply(setlist, function(i){
         !is.null(names(i));
      });
      if (all(namedCheck)) {
         setlistSigned <- setlist;
         setlist <- lapply(setlistSigned, names);
         doSigned <- TRUE;
      }
   }
   
   ## Clean up of sample sets to minimize formatting issues
   if (cleanup %in% c(TRUE)) {
      ## Set all characters to upper case
      ## Remove leading and trailing spaces
      setlist <- lapply(setlist, function(x) {
         gsub("^[ \t]+|[ \t]+$", "", toupper(x));
      });
   }
   
   ## Append object counter to retain duplicates
   if (keepdups == TRUE) {
      if (verbose) {
         printDebug("Keeping duplicates if relevant.");
      }
      dupCount <- function(setlist=setlist) {
         #count <- table(setlist);
         #paste(rep(names(count), count), unlist(sapply(count, function(x) seq(1, x))), sep=".");
         jamba::makeNames(setlist,
            suffix=".",
            renameOnes=TRUE);
      }
      mynames <- names(setlist);
      setlist <- lapply(setlist, function(x) {
         dupCount(x);
      });
      names(setlist) <- mynames;
   }
   if (doSigned) {
      setlistSigned <- lapply(jamba::nameVector(names(setlist)), function(i){
         jamba::nameVector(setlistSigned[[i]],
            setlist[[i]],
            makeNamesFunc=c);
      });
   }
   ## Create intersect matrix (removes duplicates!)
   if (verbose) {
      jamba::printDebug("overLapper(): ",
         "Creating intersect matrix.");
   }
   setunion <- sort(unique(unlist(setlist)));
   if (length(setunion) > 0) {
      if (doSigned) {
         setmatrix <- do.call(cbind, lapply(jamba::nameVector(names(setlist)), function(x) {
            #x1 <- jamba::nameVector(setunion %in% unique(setlist[[x]]),
            #   setunion);
            x1 <- jamba::nameVector(sign(setlistSigned[[x]][setunion]),
               setunion);
            x1[is.na(x1)] <- 0;
            x1;
         }));
      } else {
         setmatrix <- do.call(cbind, lapply(jamba::nameVector(names(setlist)), function(x) {
            jamba::nameVector(setunion %in% unique(setlist[[x]]),
               setunion);
         })) * 1;
      }
   } else {
      setmatrix <- matrix(ncol=length(setlist),
         nrow=0,
         dimnames=list(character(0), names(setlist)));
   }
   storage.mode(setmatrix) <- "numeric";
   
   ## Create all possible sample combinations within requested complexity levels
   if (verbose) {
      jamba::printDebug("overLapper(): ",
         "Creating pairwise comparisons.");
   }
   labels <- names(setlist);
   allcombl <- unlist(recursive=FALSE,
      lapply(complexity, function(ic) {
         combn(labels,
            m=ic,
            simplify=FALSE);
      })
   );
   names(allcombl) <- sapply(allcombl, paste, collapse=sep);
   complevels <- lengths(allcombl);
   
   ## Return intersect list for generated sample combinations
   if (type=="intersects") {
      if (verbose) {
         jamba::printDebug("overLapper(): ",
            "Handling intersects type.");
      }
      OLlist <- sapply(allcombl, function(x) {
         names(which(rowSums(abs(setmatrix[,x,drop=FALSE])) == length(x)));
      });
      if (doSigned) {
         OLfactor <- jamba::pasteByRow(setmatrix, sep=":");
      } else {
         OLfactor <- jamba::pasteByRow(setmatrix);
      }
      attr(OLfactor, "colnames") <- colnames(setmatrix);
      return(list(Set_List=setlist,
         Intersect_Matrix=setmatrix,
         Complexity_Levels=complevels,
         Intersect_List=OLlist,
         Overlap_Factor=OLfactor));
   }
   
   ## Return Venn intersect list for generated sample combinations
   if (type=="vennsets") {
      if (verbose) {
         jamba::printDebug("overLapper(): ",
            "Handling vennsets type.");
      }
      vennSets <- function(setmatrix, allcombl, index=1) {
         mycol1 <- index;
         mycol2 <- setdiff(colnames(setmatrix),index);
         cond1 <- names(which(rowSums(abs(setmatrix[,mycol1,drop=FALSE])) == length(mycol1)));
         cond2 <- names(which(rowSums(abs(setmatrix[,mycol2,drop=FALSE])) == 0));
         return(intersect(cond1, cond2));
      }
      if (verbose) {
         jamba::printDebug("overLapper(): ",
            "Handling vennsets vennOLlist.");
      }
      vennOLlist <- lapply(allcombl, function(x) {
         if (verbose) {
            jamba::printDebug(x);
         }
         vennSets(setmatrix=setmatrix,
            allcombl=allcombl,
            index=x);
      });
      if (verbose) {
         jamba::printDebug("overLapper(): ",
            "Handling vennsets vennOLfactor.");
         jamba::printDebug(jamba::cPaste(dim(setmatrix)));
      }
      if (doSigned) {
         OLfactor <- jamba::pasteByRow(setmatrix, collapse=":");
      } else {
         OLfactor <- jamba::pasteByRow(setmatrix);
      }
      attr(OLfactor, "colnames") <- colnames(setmatrix);
      return(list(Set_List=setlist,
         Intersect_Matrix=setmatrix,
         Complexity_Levels=complevels,
         Venn_List=vennOLlist,
         Overlap_Factor=OLfactor));
   }
}

#' Optimized list to signed incidence matrix
#' 
#' Optimized list to signed incidence matrix
#' 
#' This function converts a list of named vectors into
#' an incidence matrix with value for each entry
#' (row) present in each input list (column). The rows
#' are defined by the vector names, and values are
#' defined by the vector values.
#' 
#' Note that this function will store zero `0` when the input
#' vector value is zero. When this is not the desired behavior,
#' the argument `empty` can be used to distinguish missing data
#' from data that is zero, for example by setting `empty=NA`.
#' In this way a value of zero `0` indicates "present but zero",
#' and a value `NA` indicates "not present at all". This
#' distinction is helpful when comparing entities which are not
#' tested in each scenario. For example if "geneA" is present
#' and the value is `1` in one list; "geneA" is not tested in
#' the second list; therefore the absence of "geneA" of a non-zero
#' value in the second list is not counted as "non-overlapping"
#' because it was not possible for it to have a non-zero value.
#' 
#' @param setlist `list` of vectors
#' @param empty default single value used for empty/missing entries,
#'    the default `empty=0` uses zero for entries not present.
#'    Another alternative is `NA`.
#' @param do_sparse `logical` indicating whether to coerce the output
#'    to sparse matrix class `"ngCMatrix"` from the Matrix package.
#' @param coerce_sign `logical` indicating whether to coerce numeric
#'    vector values to the sign. When `coerce_sign=FALSE` the vector
#'    values are stored directly. When `coerce_sign=TRUE` the signs of
#'    the vector values are stored.
#' @param ... additional arguments are ignored.
#' 
#' @export
list2im_signed <- function
(setlist,
 empty=0,
 do_sparse=TRUE,
 force_sign=FALSE,
   ...)
{
   setnames <- lapply(setlist, names);
   setnamesunion <- Reduce("union", setnames);
   if (length(empty) == 0) {
      empty <- NA;
   } else {
      empty <- head(empty, 1);
   }
   setlistim <- do.call(cbind, lapply(setlist, function(i){
      i_match <- match(names(i), setnamesunion);
      j <- rep(empty, length(setnamesunion));
      if (force_sign) {
         j[i_match] <- sign(i);
      } else {
         j[i_match] <- i;
      }
      j;
   }))
   rownames(setlistim) <- setnamesunion;
   if (!is.character(setlistim[1,1]) &&
         do_sparse &&
         suppressPackageStartupMessages(require(Matrix))) {
      setlistim <- as(setlistim, "dgCMatrix");
   }
   return(setlistim);
}

#' Optimized list to signed incidence matrix
#' 
#' Optimized list to signed incidence matrix
#' 
#' This function rapidly converts a list of vectors into
#' an incidence matrix with value of `1` for each entry
#' (row) present in each input list (column).
#' 
#' Note that the rows in the output matrix are not sorted,
#' since this step can take several seconds when working with
#' a list whose vectors contain millions of rows.
#' 
#' @param setlist `list` of vectors
#' @param empty default single value used for empty/missing entries,
#'    the default `empty=0` uses zero for entries not present.
#'    Another alternative is `NA`.
#' @param do_sparse `logical` indicating whether to coerce the output
#'    to sparse matrix class `"ngCMatrix"` from the Matrix package.
#' @param ... additional arguments are ignored.
#' 
#' @export
list2im_opt <- function
(setlist,
 empty=0,
 do_sparse=TRUE,
 ...)
{
   setnamesunion <- Reduce("union", setlist);
   if (length(empty) == 0) {
      empty <- NA;
   } else {
      empty <- head(empty, 1);
   }
   setlistim <- do.call(cbind, lapply(setlist, function(i){
      i_match <- match(i, setnamesunion);
      j <- rep(empty,
         length(setnamesunion));
      j[i_match] <- 1;
      j;
   }))
   rownames(setlistim) <- setnamesunion;
   if (do_sparse && suppressPackageStartupMessages(require(Matrix))) {
      setlistim <- as(setlistim, "ngCMatrix");
   }
   return(setlistim);
}

#' match matrix to matrix
#' 
#' @export
match.matrix <- function
(x,
 table,
 MARGIN=1L,
 ...)
{
   temp_x <- asplit(x, MARGIN);
   temp_table <- asplit(table, MARGIN);
   match(temp_x, temp_table);
}

#' match incidence matrix to incidence matrix
#' 
#' match incidence matrix to incidence matrix
#' 
#' This function is an optimized case of `match.matrix()`
#' that requires the matrix values to be `c(-1, 0, 1)`
#' such that each row can be condensed to a single integer
#' value, and a signed variant of that integer value.
#' 
#' Once converted, the `match()` operation is quite fast.
#' 
#' There are alternatives to `match.matrix()`, first is
#' to convert the `matrix` to `data.frame`, which itself is
#' somewhat slow. Next, run `match.data.frame()` which is
#' also slow, because it ultimately runs `match()` on two
#' `list` objects.
#' 
#' @examples
#' ny <- 1000000;
#' nx <- 4;
#' x <- matrix(ncol=nx, nrow=ny, data=sample(c(-1,0,1), size=nx*ny, replace=TRUE));
#' x_u <- unique(x);
#' system.time(u_match <- match.im(x, x_u));
#' ## 1 second
#' system.time(u_match1 <- match.matrix(x, x_u));
#' ## 7 seconds
#' table(u_match1 == u_match);
#' 
#' x_df <- data.frame(x);
#' x_u_df <- data.frame(x_u);
#' system.time(u_match2 <- match(asplit(x_df, 1), asplit(x_u_df, 1)));
#' ## 10 seconds
#' table(u_match2 == u_match);
#' 
#' @export
match.im <- function
(x,
 table,
 MARGIN=1L,
 ...)
{
   ## Todo: if ncol(x) > 16 then split into chunks of 16 columns
   ## so the summary value is still an integer
   xseq <- 2^(seq_len(ncol(x))-1);
   x_score <- colSums(t(abs(x)) * xseq);
   x_score_signed <- colSums(t(x) * xseq);
   x_key <- paste(x_score, x_score_signed);

   y_score <- colSums(t(abs(table)) * xseq);
   y_score_signed <- colSums(t(table) * xseq);
   y_key <- paste(y_score, y_score_signed);
   
   match(x_key, y_key);
}

