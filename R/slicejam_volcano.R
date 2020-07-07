
## volcano plot functions

#' Volcano plot
#' 
#' Volcano plot
#' 
#' @export
volcanoPlot <- function
(x,
   n=NULL,
   fcCol=head(jamba::vigrep(c("log.*fc|fold.*ch","fc","fold"), colnames(x)),1),
   pvCol=head(jamba::provigrep(c("^padj|adj.*p","fdr","q.*val","p.*val"), colnames(x)),1),
   intensityCols=jamba::provigrep(c("maxgroupmean|maxmean","ave.*expr","^fpkm"), colnames(x)),
   pvalueCutoff=0.01,
   fcCutoff=0,
   fcCutoffInLog=FALSE,
   intensityCutoff=NULL,
   main="Volcano Plot",
   lineMain=3,
   subMain=NULL,
   font.main=1,
   symmetricAxes=TRUE,
   ptCex=0.9,
   ptPch=21,
   xlab=paste("Fold change (", fcCol, ")", sep=""),
   ylab=paste("Significance (", pvCol, ")", sep=""),
   doStatsSubtitles=TRUE, statsLine=4, statsCex=0.8,
   minPvalue=1e-20,
   pvalueFloor=minPvalue,
   maxLog2FC=10,
   nXlabels=13,
   nYlabels=7,
   xRange=NULL,
   yRange=NULL,
   fcCeiling=NULL,
   hitType="hits",
   hitColorSet=c("#77777777", "#99000055"),
   pointColors=NULL,
   pointBgColors=NULL,
   pointColorSet=c("base"="#77777733",
      "up"="#99000055",
      "down"="#00009955",
      #"base.highlight"="#EEEE55DD",
      "base.highlight"="#FFDD55DD",
      "up.highlight"="#FFDD55DD",
      "down.highlight"="#FFDD55DD"),
   pointBorderSet=c("base"="#33333333",
      "up"="#77000088",
      "down"="#00007788",
      "base.highlight"="#773322DD",
      "up.highlight"="#773322DD",
      "down.highlight"="#773322DD"),
   contrastFactor=1.5,
   usePointColorSet=TRUE,
   ablineColor="#000000AA",
   highlightPoints=NULL,
   highlightHits=FALSE,
   highlightItemType=paste("highlighted", hitType, sep=" "),
   highlightCex=1,
   pointsToLabel=NULL, noObscureMode=TRUE, initialAngle=0, initialRadius=0,
   labelCex=1, labelMar=labelCex*1.1, boxColor="#EEBB7788", boxBorderColor="#000000AA",
   labelFixedCoords=NULL,
   doHighlightLabels=TRUE,
   sectionLabelSpacing=12, labelN=30,
   doJitter=FALSE,
   jitterAmount=NULL,
   geneColumn=vigrep("gene|symbol", colnames(x))[1], #NULL,
   hitsByGene=FALSE,
   printMultiGeneHits=FALSE,
   doSmoothScatter=TRUE,
   smoothScatterFunc=plotSmoothScatter,
   #smoothColramp=colorRampPalette(c("white", "lightblue", "blue", "orange", "orangered2")),
   smoothColramp=colorRampPalette(c("white", "lightblue", "lightskyblue3", "royalblue", "darkblue", "orange", "darkorange1", "orangered2")),
   #smoothColramp=colorRampPalette(c("white", "blue", "orange", "red")),
   useRaster=TRUE,
   doBoth=FALSE, 
   bothColorOnlyHits=TRUE,
   orientation="up",
   addPlot=FALSE, 
   pointSubset=NULL, 
   plotOnlySubset=FALSE,
   transformation=function(x){x^0.24;}, # use 0.14 for very large datasets
   nbin=256, nrpoints=0,
   doBlockArrows=TRUE, blockArrowFont=1,
   blockArrowHitColor="#E67739FF", blockArrowUpColor="#990000FF", blockArrowDownColor="#000099FF",
   blockArrowLabelHitColor="#FFFFAA", blockArrowLabelUpColor="#FFFFAA", blockArrowLabelDownColor="#FFFFAA",
   blockArrowCex=c(1,1), blockArrowLabelCex=1,
   doShadowTextArrows=TRUE,
   doTopHist=FALSE, doTopHistCutoffs=c("pvalue", "foldchange"),
   topHistBreaks=100, topHistColor="#000099FF",
   topHistCexSub=1, topHistSubLine=0,
   topHistPlotFraction=1/5, topHistBy=0.20,
   cexSub=0.8, cexMain=1.5, cex.axis=1.2,
   lineXlab=2.5, lineYlab=4,
   labelHits=!doBlockArrows, verbose=FALSE,
   ...)
{

   ## Purpose is to wrapper a simple volcano plot with log-friendly axis labels
   ##
   ## orientation can be "up" for upright, "right" to tilt 90 degree on its side
   ##
   ## minPvalue is designed to set NA or '0' values to a fixed value
   ## pvalueFloor is designed to set extremely low P-values (1e-150) to a minimum
   ## sectionLabelSpacing is the fraction of the y-axis range away from each border
   ##    used to position the highlight hit counts labels (if highlightPoints is not NULL)
   ##
   ## A change was made to use pointColorSet to set the colors, instead of using
   ## hitColorSet, so the colors could independently be manipulated.
   ## To use this scheme, set usePointColorSet=TRUE
   ##
   ## The graph plots fold change on the x-axis, but with log2 scale, the numbers
   ## are technically being reported in normal space. However, to give it a log2 axis label:
   ## xlab=expression(log[2]*"(Fold change)"
   blockArrowCex <- rep(blockArrowCex, length.out=2);
   blockArrowCex <- rep(blockArrowCex, length.out=2);
   if (!usePointColorSet) {
      pointColorSet["up"] <- hitColorSet[2];
      pointColorSet["down"] <- hitColorSet[2];
      pointColorSet["base"] <- hitColorSet[1];
      pointColorSet["base.highlight"] <- hitColorSet[2];
      pointColorSet["up.highlight"] <- hitColorSet[2];
      pointColorSet["down.highlight"] <- hitColorSet[2];
      if (is.null(pointBgColors)) {
         pointAlphas <- (jamba::col2alpha(pointColorSet)*2+1)/3;
         if (verbose) {
            jamba::printDebug("pointAlphas: ",
               format(digits=2, pointAlphas),
               c("orange", "seagreen"));
         }
         pointDFs <- ifelse(jamba::col2hcl(pointColorSet)["L",] <= 60, -1.5, 1.5);
         pointBorderSet <- jamba::makeColorDarker(pointColorSet,
            fixAlpha=pointAlphas,
            darkFactor=pointDFs,
            sFactor=1.5);
      }
   } else {
      ## Only alter parameters which were provided in the function call,
      ## but keep default values which were not defined
      pointColorSet <- update_function_params(functionName="volcanoPlot",
         param_name="pointColorSet",
         new_values=pointColorSet);
      pointBorderSet <- update_function_params(functionName="volcanoPlot",
         param_name="pointBorderSet",
         new_values=pointBorderSet);
      if (verbose) {
         jamba::printDebug(" pointColorSet:");
         jamba::printDebugI(pointColorSet);
         jamba::printDebug("pointBorderSet:");
         jamba::printDebugI(pointBorderSet);
      }
   }
   ##
   if (doBoth && bothColorOnlyHits) {
      hitColorSet[1] <- "#FFFFFF00";
      pointColorSet["base"] <- hitColorSet[1];
   }
   
   if (!pvCol %in% colnames(x)) {
      pvCol <- head(jamba::provigrep(gsub("[-()_.]+", ".+", pvCol), colnames(x)), 1);
   }
   if (!fcCol %in% colnames(x)) {
      fcCol <- head(jamba::provigrep(gsub("[-()_.]+", ".+", fcCol), colnames(x)), 1);
   }
   if (verbose) {
      jamba::printDebug("volcanoPlot(): ",
         "fcCol:",
         fcCol);
      jamba::printDebug("volcanoPlot(): ",
         "pvCol:",
         pvCol);
   }
   if (length(pvCol) == 0 || length(fcCol) == 0) {
      stop("volcanoPlot(): pvCol and fcCol must match colnames(x).");
   }
   origPar <- par();
   
   ## Allow processing only a subset 'n' rows of data
   if (is.null(n)) {
      n <- nrow(x);
   } else if (n < nrow(x)) {
      x <- x[1:n,,drop=FALSE];
   }
   
   if (is.null(rownames(x)) || length(jamba::tcount(rownames(x), minCount=2)) > 0) {
      rownames(x) <- jamba::makeNames(rep("row", n));
   }

   if (any(intensityCols %in% colnames(x))) {
      intensityCols <- head(intersect(intensityCols, colnames(x)), 1);
   } else {
      intensityCols <- head(jamba::provigrep(intensityCols,
         colnames(x)), 1);
   }
   if (length(intensityCutoff) > 0) {
      if (verbose) {
         printDebug("volcanoPlot(): ",
            "intensityCols:",
            intensityCols);
      }
      if (length(intensityCols) > 0) {
         if (verbose) {
            printDebug("Applying intensityCutoff.");
         }
         if (!suppressPackageStartupMessages(require(matrixStats, quietly=TRUE))) {
            stop("The matrixStats package is required.");
         }
         intensityValues <- rowMaxs(as.matrix(x[,intensityCols,drop=FALSE]));
         metIntensity <- intensityValues > intensityCutoff;
         if (verbose) {
            print(table(metIntensity));
         }
      } else {
         metIntensity <- rep(TRUE, nrow(x));
      }
   } else {
      metIntensity <- rep(TRUE, nrow(x));
   }
   pvHitsUpWhich <- which(x[,pvCol] <= pvalueCutoff & (x[,fcCol]) >= fcCutoff & metIntensity);
   pvHitsDownWhich <- which(x[,pvCol] <= pvalueCutoff & -(x[,fcCol]) >= fcCutoff & metIntensity);
   pvHitsBothWhich <- which(x[,pvCol] <= pvalueCutoff & abs(x[,fcCol]) >= fcCutoff & metIntensity);
   pvHitsUp <- length(pvHitsUpWhich);
   pvHitsDown <- length(pvHitsDownWhich);
   pvHits <- length(pvHitsBothWhich);
   upHits <- rownames(x)[pvHitsUpWhich];
   downHits <- rownames(x)[pvHitsDownWhich];

   ## Define point colors
   if (is.null(pointColors)) {
      pvColors1 <- jamba::nameVector(rep(pointColorSet["base"], n), rownames(x));
      pvColors1[upHits] <- pointColorSet["up"];
      pvColors1[downHits] <- pointColorSet["down"];
   } else {
      pvColors1 <- pointColors;
   }
   if (is.null(pointBgColors)) {
      #pvBgColors1 <- makeColorDarker(pvColors1, darkFactor=1.5);
      pvBgColors1 <- nameVector(rep(pointBorderSet["base"], n), rownames(x));
      pvBgColors1[upHits] <- pointBorderSet["up"];
      pvBgColors1[downHits] <- pointBorderSet["down"];
   } else {
      pvBgColors1 <- pointBgColors;
   }
   if (verbose) {
      printDebug("head(unique(pvColors1)): ");
      printDebug(head(unique(pvColors1)));
   }
   if (class(x[,pvCol]) %in% c("numeric", "integer")) {
      pvValues <- nameVector(x[,pvCol], rownames(x)[1:n]);
   } else {
      pvValues <- nameVector(as.numeric(x[,pvCol]), rownames(x)[1:n]);
   }
   pvValues[is.na(pvValues)] <- 1;
   if (verbose) {
      printDebug("minPvalue: ", minPvalue);
      printDebug("pvValues: ", head(pvValues));
      printDebug("any(pvValues < minPvalue): ", any(pvValues < minPvalue));
   }
   if (!is.null(minPvalue) && any(pvValues < minPvalue)) {
      pvValues[pvValues < minPvalue] <- minPvalue;
   }
   if (!is.null(pvalueFloor) && any(pvValues < pvalueFloor)) {
      pvValues[pvValues < pvalueFloor] <- pvalueFloor;
   }
   yValues1 <- -log10(pvValues);
   
   fcValues <- nameVector(x[,fcCol], rownames(x));
   
   fcValues[is.na(fcValues)] <- 0;
   fcValues[abs(fcValues) > maxLog2FC] <- sign(fcValues[abs(fcValues) > maxLog2FC]) * maxLog2FC;
   xValues1 <- fcValues;

   
   ## Allow highlighting a subset of points, which overrides the pointSubset argument
   if (is.null(highlightPoints) && highlightHits) {
      if (verbose) {
         printDebug("   volcanoPlot defining highlightPoints as up/down hits.");
      }
      highlightPoints <- rownames(x)[unique(c(pvHitsUpWhich, pvHitsDownWhich))];
      highlightItemType <- hitType;
   }
   if (!is.null(pointsToLabel)) {
      row2name <- c(jamba::nameVector(rownames(x), x[,geneColumn]),
         jamba::nameVector(rownames(x)));
      pointsToLabel <- jamba::flipVector(row2name[pointsToLabel],
         makeNamesFunc=c);
   }
   if (verbose) {
      jamba::printDebug("   volcanoPlot length(highlightPoints):", 
         jamba::formatInt(length(highlightPoints)));
   }
   if (verbose) {
      jamba::printDebug("   doBoth=", doBoth);
   }
   if (!is.null(highlightPoints)) {
      if (verbose) {
         jamba::printDebug("Creating subset of points for highlighting.");
      }
      if (is.null(doSmoothScatter) || doSmoothScatter) {
         if (verbose) {
            jamba::printDebug("Setting doBoth=", "TRUE");
         }
         doBoth <- TRUE;
      }
      pointSubset <- which(rownames(x) %in% highlightPoints | x[,geneColumn] %in% highlightPoints);
      if (verbose) {
         jamba::printDebug("   volcanoPlot head(pointSubset):", head(pointSubset));
         jamba::printDebug("   volcanoPlot length(pointSubset):",
            jamba::formatInt(length(pointSubset)));
      }
      ## If we're plotting all individual points, then we can color-code the hits, and the highlighted
      ## items all independently
      highlightUpHits <- intersect(rownames(x)[pointSubset], upHits);
      highlightDownHits <- intersect(rownames(x)[pointSubset], downHits);
      highlightBaseHits <- setdiff(rownames(x)[pointSubset], c(upHits, downHits));
      pvColors1[highlightBaseHits] <- pointColorSet["base.highlight"];
      pvColors1[highlightUpHits] <- pointColorSet["up.highlight"];
      pvColors1[highlightDownHits] <- pointColorSet["down.highlight"];
      pvBgColors1[highlightBaseHits] <- pointBorderSet["base.highlight"];
      pvBgColors1[highlightUpHits] <- pointBorderSet["up.highlight"];
      pvBgColors1[highlightDownHits] <- pointBorderSet["down.highlight"];
      
      pointSubsetLabels <- x[pointSubset,geneColumn];
      if (verbose) {
         printDebug("highlightUpHits: ", head(highlightUpHits), c("orange", "lightblue"));
         printDebug("highlightDownHits: ", head(highlightDownHits), c("orange", "lightblue"));
         printDebug("highlightBaseHits: ", head(highlightBaseHits), c("orange", "lightblue"));
         printDebug("pointSubsetLabels: ", head(pointSubsetLabels), c("orange", "lightblue"));
         printDebug("pvColors1[highlightDownHits]:");
         printDebug(names(pvColors1[highlightDownHits]), list(pvColors1[highlightDownHits]));
         printDebug("pvBgColors1[highlightDownHits]:");
         printDebug(names(pvBgColors1[highlightDownHits]), list(pvBgColors1[highlightDownHits]));
      }
      names(pointSubset) <- names(xValues1)[pointSubset];
      if (length(highlightCex) == 1) {
         pointSubsetCex <- highlightCex;
      } else {
         ## Order the given highlightCex alongside the matching rows in the data
         names(highlightCex) <- highlightPoints;
         pointSubsetCex <- highlightCex[names(pointSubset)];
      }
      if (length(pointSubset) == 0) {
         pointSubset <- NULL;
         pointSubsetCex <- 0;
         if (verbose) {
            printDebug("No points found to highlight.", fgText="yellow");
         }
      } else {
         ## Generate hit counts among the highlighted points
         pvHighlightHitsUp <- length(which(x[1:n,pvCol][pointSubset] <= pvalueCutoff & (x[1:n,fcCol][pointSubset]) >= fcCutoff));
         pvHighlightHitsDown <- length(which(x[1:n,pvCol][pointSubset] <= pvalueCutoff & -(x[1:n,fcCol][pointSubset]) >= fcCutoff));
         pvHighlightHitsOther <- length(pointSubset) - (pvHighlightHitsUp + pvHighlightHitsDown);
         if (verbose) {
            printDebug("   pvHighlightHitsUp: ", pvHighlightHitsUp, c("orange", "lightblue"));
            printDebug("   pvHighlightHitsDown: ", pvHighlightHitsDown, c("orange", "lightblue"));
            printDebug("   pvHighlightHitsOther: ", pvHighlightHitsOther, c("orange", "lightblue"));
         }
      }
   }
   if (!is.null(pointSubset) && plotOnlySubset) {
      xValues1 <- xValues1[pointSubset];
      yValues1 <- yValues1[pointSubset];
      pvColors1 <- pvColors1[pointSubset];
      pvBgColors1 <- pvBgColors1[pointSubset];
      pointSubset <- NULL;
      pointSubsetCex <- 0;
   }

   
   ## Switch pvBgColors1 and pvColors1
   pvColors2 <- pvColors1;
   pvColors1 <- pvBgColors1;
   pvBgColors1 <- pvColors2;
   
   if (!is.null(fcCeiling)) {
      xValues1[abs(xValues1) > fcCeiling] <- sign(xValues1[abs(xValues1) > fcCeiling]) * fcCeiling;
   }
   
   if (doJitter) {
      if (verbose) {
         printDebug("Jittering points by:", format(digits=2, jitterAmount));
      }
      xValues2 <- xValues1;
      yValues2 <- yValues1;
      xValues1 <- jitter2(xValues1, factor=as.numeric(doJitter), amount=jitterAmount);
      yValues1 <- jitter2(yValues1, factor=as.numeric(doJitter), amount=jitterAmount);
   }
   
   if (is.null(xRange)) {
      xRange <- range(xValues1);
   }
   if (is.null(yRange)) {
      yRange <- range(c(0, yValues1, -log10(pvalueCutoff/5)));
   }
   if (symmetricAxes) {
      xRange <- c(-max(abs(xRange)), max(abs(xRange)));
   }
   
   if (!is.null(pointSubset)) {
      ###############################################
      ## Position labels in quadrants of the plot
      namesX <- c(xRange[1] - xRange[1]/sectionLabelSpacing, xRange[2] - xRange[2]/sectionLabelSpacing);
      namesY <- rep(yRange[2] - yRange[2]/sectionLabelSpacing, 2);
      namesLabel <- paste(c(format(big.mark=",", pvHighlightHitsDown), format(big.mark=",", pvHighlightHitsUp)), highlightItemType, sep=" ");
      namesLabel <- gsub("^(1 .+)s$", "\\1", namesLabel);
      if (verbose) {
         printDebug("namesLabel: ", namesLabel, c("orange", "seagreen"));
      }
   }
   
   parList <- list();
   parList[["prePlots"]] <- par();

   
   ## labelCoords will have the return data from addNonOverlappingLabels() but only
   ## if we end up calling that method
   labelCoords <- NULL;
   ##########################################################
   ## Histogram along the top border, set up the details here
   if (doTopHist) {
      if (length(grep("pval", doTopHistCutoffs)) > 0) {
         ## Apply P-value filtering
         pvHitsWhich <- which(x[1:n,pvCol] <= pvalueCutoff);
      } else {
         pvHitsWhich <- 1:n;
      }
      if (length(grep("fc|fold", doTopHistCutoffs)) > 0) {
         ## Apply fold change filtering
         fcHitsWhich <- which(abs(x[1:n,fcCol]) >= fcCutoff);
      } else {
         fcHitsWhich <- 1:n;
      }
      histWhich <- intersect(pvHitsWhich, fcHitsWhich);
      #return(histWhich);
      
      ## Color the bars consistent with the block arrows
      topBarColors <- c("#000099FF", "#990000FF");
      parMar <- par("mar");
      topHistBreaks <- seq(from=xRange[1], to=xRange[2], by=topHistBy);
      if (max(topHistBreaks) < xRange[2]) {
         if (verbose) {
            printDebug("Fixing topHistBreaks: ", cPaste(topHistBreaks), c("orange", "lightblue"));
         }
         topHistBreaks <- c(topHistBreaks, xRange[2]);
      }
      topHist <- hist(xValues2[histWhich], breaks=topHistBreaks, plot=FALSE);
      plotZones <- matrix(c(1,2), nrow=2);
      layout(plotZones, widths=c(1), heights=c(topHistPlotFraction, 1-topHistPlotFraction));
      ## Change margins, then plot the top histogram
      ## default is par(mar=c(5.1, 4.1, 4.1, 2.1));
      par("mar"=c(1, parMar[2], parMar[3], parMar[4]));
      parList[["preTopHist"]] <- par();
      if (verbose) {
         printDebug("par(mgp): ", cPaste(par("mgp")), c("orange", "lightblue"));
         printDebug("length(topHistBreaks): ", length(topHistBreaks), c("orange", "lightblue"));
      }
      r1 <- as.integer(length(topHistBreaks)/2);
      topHistCol <- rep(topBarColors, c(r1, r1+1));
      barplot(topHist$counts, axes=FALSE, ylim=c(0, max(topHist$counts)), space=0,
         col=topHistCol, border=makeColorDarker(topHistCol), horiz=FALSE, las=2, cex.axis=1);
      prettyAt1 <- pretty(c(0, max(topHist$counts)), n=4);
      prettyAt1 <- prettyAt1[prettyAt1 <= max(topHist$counts)];
      axis(2, las=2, cex.axis=1.3, at=prettyAt1, ...);
      title(xlab=paste("Distribution of ", hitType, sep=""), cex.lab=topHistCexSub, xpd=TRUE, outer=FALSE, line=topHistSubLine);
      parList[["postTopHist"]] <- par();
      par("mar"=c(parMar[1], parMar[2], parMar[3]-1.5, parMar[4]));
   }

   
   ######################
   ## Smooth scatter plot
   if (doSmoothScatter || doBoth) {
      ## Upright volcano plot (standard orientation)
      smoothScatterFunc(x=xValues1, y=yValues1, #col=pvColors1, bg=pvBgColors1,
         xaxt="n", yaxt="n", xlab="", ylab="", transformation=transformation, useRaster=useRaster,
         xlim=xRange, ylim=yRange, nbin=nbin, colramp=smoothColramp, add=addPlot, nrpoints=nrpoints, ...);
      title(xlab=xlab, line=lineXlab, cex.axis=cex.axis, ...);
      title(ylab=ylab, line=lineYlab, cex.axis=cex.axis, ...);
      parUsr <- par("usr");
      if (verbose) {
         printDebug("xlab: ", xlab, c("orange", "lightblue"));
      }
      if (doBoth) {
         if (!is.null(pointSubset)) {
            if (verbose) {
               printDebug("Following doBoth to label pointSubset.");
            }
            points(x=xValues1[pointSubset], y=yValues1[pointSubset], pch=ptPch, cex=pointSubsetCex,#ptCex,
               col=pvColors1[pointSubset], bg=pvBgColors1[pointSubset], xaxt="n", yaxt="n",
               xlim=xRange, ylim=yRange, ...);
            ## Optionally labels the highlighted points, using the addNonOverlappingLabels() function
            if (doHighlightLabels) {
               printDebug("length(pointSubset):", length(pointSubset), ", labelN:", labelN);
               labelCoords <- addNonOverlappingLabels(x=xValues1[pointSubset], y=yValues1[pointSubset], initialAngle=initialAngle,
                  txt=pointSubsetLabels, n=labelN,
                  labelCex=labelCex,  labelMar=labelMar, boxColor=boxColor,
                  boxBorderColor=boxBorderColor, initialRadius=initialRadius,
                  fixedCoords=labelFixedCoords,
                  ...);
            }
            ## Add text labels indicating the number of highlighted points
            textAdjX <- c(0, 1);
            if (labelHits) {
               for(i in seq_along(namesLabel)) {
                  printDebug("   namesLabel[i]: ", namesLabel[i], c("orange", "seagreen"));
                  text(x=namesX[i], y=namesY[i], label=namesLabel[i], adj=c(textAdjX[i], 0.5));
               }
            }
         } else {
            if (verbose) {
               printDebug("Following doBoth but pointSubset is NULL.");
            }
            points(x=xValues1, y=yValues1, pch=ptPch, cex=ptCex,
               col=pvColors1, bg=pvBgColors1, xaxt="n", yaxt="n",
               xlim=xRange, ylim=yRange, ...);
         }
      }
      parList[["postSmoothScatter"]] <- par();
   }
   
   
   ## Non-smooth scatter points
   if (!doSmoothScatter || doBoth) {
      ## Non-smoothScatter
      ## Add points to an existing plot
      if (addPlot || doBoth) {
         if (!is.null(pointSubset) && !plotOnlySubset) {
            points(y=yValues1[pointSubset],
               x=xValues1[pointSubset], 
               pch=ptPch, 
               cex=pointSubsetCex,#ptCex,
               col=pvColors1[pointSubset], 
               bg=pvBgColors1[pointSubset], 
               xaxt="n", 
               yaxt="n",
               xlim=xRange, 
               ylim=yRange, 
               ...);
            if (verbose) {
               jamba::printDebug("pvColors1[pointSubset]:");
               jamba::printDebug((pvColors1[pointSubset]), list(pvBgColors1[pointSubset]));
               jamba::printDebug("pvBgColors1[pointSubset]:");
               jamba::printDebug((pvBgColors1[pointSubset]), list(pvBgColors1[pointSubset]));
            }
         } else {
            points(x=xValues1, 
               y=yValues1, 
               pch=ptPch, 
               cex=ptCex,
               col=pvColors1, 
               bg=pvBgColors1, 
               xaxt="n",
               yaxt="n",
               xlim=xRange, 
               ylim=yRange, 
               ...);
         }
      } else {
         ## Create new all-point scatter plot
         doPoints <- which(!names(xValues1) %in% pointSubset);
         if (verbose) {
            printDebug("length(doPoints): ", length(doPoints), c("orange", "lightblue"));
         }
         plot(x=xValues1[doPoints], y=yValues1[doPoints], pch=ptPch, cex=ptCex,
            col=pvColors1[doPoints], bg=pvBgColors1[doPoints], xaxt="n", yaxt="n",
            xlab="", ylab="",
            xlim=xRange, ylim=yRange, ...);
         title(xlab=xlab, line=lineXlab, cex.axis=cex.axis, ...);
         title(ylab=ylab, line=lineYlab, cex.axis=cex.axis, ...);
         if (verbose) {
            printDebug("ping", "orange");
         }
         ######################
         ## Optionally label a specific subset of points
         if (!is.null(pointsToLabel) && length(pointsToLabel) > 0) {
            if (verbose) {
               printDebug("Labeling ", length(pointsToLabel), " pointsToLabel.", c("orange", "lightblue"));
            }
            ## If pointsToLabel has names, they're used to match rownames(x), and the values
            ## are used as text labels on the plot.  This allows for adding a label which is
            ## different from the data matrix rownames.
            if (is.null(names(pointsToLabel))) {
               pointsToLabel <- nameVector(pointsToLabel);
            }
            labelCoords <- addNonOverlappingLabels(x=xValues1[names(pointsToLabel)], y=yValues1[names(pointsToLabel)],
               labelCex=labelCex, labelMar=labelMar, boxColor=boxColor, boxBorderColor=boxBorderColor,
               initialRadius=initialRadius,
               txt=pointsToLabel, noObscureMode=noObscureMode, initialAngle=initialAngle,
               fixedCoords=labelFixedCoords);
         }
         ######################
         ## Borrowed from above
         if (!is.null(pointSubset)) {
            if (verbose) {
               printDebug("length(pointSubset): ", length(pointSubset), c("orange", "lightblue"));
            }
            points(x=xValues1[pointSubset], y=yValues1[pointSubset], pch=ptPch, cex=pointSubsetCex,#ptCex,
               col=pvColors1[pointSubset], bg=pvBgColors1[pointSubset], xaxt="n", yaxt="n",
               xlim=xRange, ylim=yRange, cex.axis=cex.axis, ...);
            ## Optionally labels the highlighted points, using the addNonOverlappingLabels() function
            if (doHighlightLabels) {
               labelCoords <- addNonOverlappingLabels(x=xValues1[pointSubset], y=yValues1[pointSubset],
                  labelCex=labelCex, labelMar=labelMar, boxColor=boxColor, boxBorderColor=boxBorderColor,
                  initialRadius=initialRadius,
                  txt=pointSubsetLabels, n=labelN, initialAngle=initialAngle,
                  fixedCoords=labelFixedCoords);
            }
            ## Add text labels indicating the number of highlighted points
            textAdjX <- c(0, 1);
            if (labelHits) {
               for(i in seq_along(namesLabel)) {
                  text(x=namesX[i], y=namesY[i], label=namesLabel[i], adj=c(textAdjX[i], 0.5));
               }
            }
         } else {
            points(x=xValues1, y=yValues1, pch=ptPch, cex=ptCex,
               col=pvColors1, bg=pvBgColors1, xaxt="n", yaxt="n",
               xlim=xRange, ylim=yRange, ...);
         }
         #######
      }
      parList[["postScatter"]] <- par();
   }
   
   
   multiGenesUp <- character(0);
   multiGenesDown <- character(0);
   if (!addPlot) {
      #labelUp <- paste(format(big.mark=",", pvHitsUp), "up-regulated", hitType);
      #labelUp <- gsub("((^|[ ])1 .+)s$", "\\1", labelUp);
      #labelDown <- paste(format(big.mark=",", pvHitsDown), "down-regulated", hitType);
      #labelDown <- gsub("((^|[ ])1 .+)s$", "\\1", labelDown);
      
      ## Define labels
      labelUp <- paste(format(big.mark=",", pvHitsUp), hitType, "up");
      labelDown <- paste(format(big.mark=",", pvHitsDown), hitType, "down");
      
      ## Get rid of "1 genes", instead a more civilized "1 gene"
      labelUp <- gsub("((^|[ ])1 .+)s$", "\\1", labelUp);
      labelDown <- gsub("((^|[ ])1 .+)s$", "\\1", labelDown);
      if (!is.null(geneColumn) && hitsByGene && geneColumn %in% colnames(x)) {
         ## Aggregate hit counts by gene instead of by assay/probe/whatever
         if (verbose) {
            printDebug("Aggregating hit totals by gene column '", geneColumn, "'", c("orange", "lightblue"));
         }
         pvHitsUp1 <- x[(which(x[1:n,pvCol] <= pvalueCutoff & (x[1:n,fcCol]) >= fcCutoff)),geneColumn];
         multiGenesUp <- tcountPaste(pvHitsUp1[pvHitsUp1 %in% names(tcount(pvHitsUp1)[tcount(pvHitsUp1) > 1])]);
         pvHitsDown1 <- x[(which(x[1:n,pvCol] <= pvalueCutoff & -(x[1:n,fcCol]) >= fcCutoff)),geneColumn];
         multiGenesDown <- tcountPaste(pvHitsDown1[pvHitsDown1 %in% names(tcount(pvHitsDown1)[tcount(pvHitsDown1) > 1])]);
         
         if (printMultiGeneHits) {
            title(adj=0, sub=paste("Genes Multiply-Up:", wordWrap(multiGenesUp, 70, cleanText=FALSE, justify="left"), "\n",
               "Genes Multiply-Down:", wordWrap(multiGenesDown, 70, cleanText=FALSE, justify="left"), sep=""),
               cex.sub=0.6, line=-2, outer=TRUE);
         }
         pvHitsUp <- length(unique(pvHitsUp1));
         pvHitsDown <- length(unique(pvHitsDown1));
         hitType <- "genes";
         labelUp1 <- labelUp;
         labelDown1 <- labelDown;
         labelUp <- paste("(", format(big.mark=",", pvHitsUp), hitType, ")");
         labelUp <- gsub("((^|[^0-9])1 .+)s$", "\\1", labelUp);
         labelUp <- paste(labelUp1, labelUp, sep="\n");
         labelDown <- paste("(", format(big.mark=",", pvHitsDown), hitType, ")");
         labelDown <- gsub("((^|[^0-9])1 .+)s$", "\\1", labelDown);
         labelDown <- paste(labelDown1, labelDown, sep="\n");
      }
      labelX <- c(range(xValues1)[1] - range(xValues1)[1]/5,
         range(xValues1)[2] - range(xValues1)[2]/5);
      labelY <- rep(range(yValues1)[2] - range(xValues1)[2]/5, 2);
      if (!addPlot && labelHits && !doBlockArrows) {
         #text(x=labelX,
         #   y=labelY,
         #   label=c(labelDown, labelUp));
         jamba::drawLabels(preset=c("topleft", "topright"),
            txt=c(labelDown, labelUp),
            labelCex=1,
            drawBox=FALSE);
      }
      
      ## Block Arrows
      parList[["preBlockArrows"]] <- par();
      par("xpd"=FALSE);
      
      if (verbose) {
         jamba::printDebug("y-axis at:", unique(as.integer(pretty(yRange, n=nYlabels))));
         jamba::printDebug("x-axis at:", unique(as.integer(pretty(xRange, n=nXlabels))));
      }
      logAxis(2, 
         at=unique(as.integer(pretty(yRange, n=nYlabels))),
         value=FALSE, 
         base=10, 
         makeNegative=TRUE, 
         cex.axis=cex.axis*0.8,
         ...);
      #logAxis(1, at=unique(as.integer(pretty(xRange, n=nXlabels))),
      #   value=TRUE, base=2, cex.axis=cex.axis*0.8, ...);
      jamba::minorLogTicksAxis(1,
         logBase=2,
         displayBase=2,
         majorCex=cex.axis,
         minorCex=cex.axis*0.7,
         symmetricZero=TRUE,
         offset=0,
         ...);
      abline(h=-log10(pvalueCutoff), 
         v=unique(c(-fcCutoff, fcCutoff)),
         lty="dashed", 
         col=ablineColor);
      if (doBlockArrows) {
         hitCol <- hsv(h=0.06, 
            s=0.75, 
            v=0.9, 
            alpha=1);
         if (doTopHist) {
            rightAdj <- 1.2;
         } else {
            rightAdj <- 1;
         }
         blockArrowMargin(axisPosition="rightAxis",
            ybottom=-log10(pvalueCutoff),
            arrowPosition="top",
            doShadowText=doShadowTextArrows,
            blockWidthPercent=5*blockArrowCex[2]*rightAdj, 
            arrowLengthPercent=3*blockArrowCex[2]*rightAdj,
            arrowLabel=paste(format(big.mark=",", pvHits), " significant ", hitType, sep=""),
            col=blockArrowHitColor,
            labelCex=blockArrowLabelCex,
            labelFont=blockArrowFont, 
            arrowLabelColor=blockArrowLabelHitColor,
            arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelHitColor), 0.3));
         blockArrowMargin(axisPosition="topAxis",
            xleft=fcCutoff, 
            arrowPosition="right",
            doShadowText=doShadowTextArrows,
            blockWidthPercent=5*blockArrowCex[1], 
            arrowLengthPercent=3*blockArrowCex[1],
            arrowLabel=labelUp, 
            col=blockArrowUpColor, 
            labelCex=blockArrowLabelCex,
            labelFont=blockArrowFont, 
            arrowLabelColor=blockArrowLabelUpColor,
            arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelUpColor), 0.3));
         blockArrowMargin(axisPosition="topAxis", 
            xright=-fcCutoff, 
            arrowPosition="left",
            doShadowText=doShadowTextArrows,
            blockWidthPercent=5*blockArrowCex[1], 
            arrowLengthPercent=3*blockArrowCex[1],
            arrowLabel=labelDown, 
            col=blockArrowDownColor, 
            labelCex=blockArrowLabelCex,
            labelFont=blockArrowFont, 
            arrowLabelColor=blockArrowLabelDownColor,
            arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelDownColor), 0.3));
      }

      parList[["postBlockArrows"]] <- par();
      
      ## Display the significance cutoff used
      #subTitle <- paste("Significance cutoff <= ", pvalueCutoff, ", and fold change cutoff > ", round(digits=2, 2^fcCutoff));
      if (!fcCutoffInLog) {
         fcCutoffLabel <- signif(digits=2, 2^fcCutoff);
      } else {
         fcCutoffLabel <- signif(digits=2, fcCutoff);
      }
      if (doStatsSubtitles) {
         if (fcCutoff > 0) {
            #subTitle <- paste(sep="", "Significance cutoff <= ", pvalueCutoff, ", and ",
            #                  ifelse(fcCutoffInLog, "log2 ", ""),
            #                  "fold change cutoff > ", fcCutoffLabel);
            subTitle <- paste(sep="", "Significance cutoff <= ", pvalueCutoff, ", and ",
               ifelse(fcCutoffInLog, "log2 ", ""),
               "fold change cutoff >= ", fcCutoffLabel);
         } else {
            subTitle <- paste("Significance cutoff <= ", pvalueCutoff);
         }
         title(sub=subTitle, cex.sub=statsCex, line=statsLine, ...);
         
         ## Display the total points
         totalPointsSub <- paste("Total points: ", format(big.mark=",", nrow(x)));
         parMar <- par("mar");
         par("mar"=c(parMar[1], 0.5, parMar[3:4]));
         title(sub=totalPointsSub, adj=0.01, cex.sub=statsCex, line=statsLine, ...);
         par("mar"=parMar);
      }
      
      ## Overall Title
      if (doTopHist && !orientation %in% c("right")) {
         #par(parList[["postTopHist"]]);
         #par("mfg"=c(1,1));
         #par("mar"=parList[["postTopHist"]]$mar);
         #par("plt"=parList[["postTopHist"]]$plt);
         #par("usr"=parList[["postTopHist"]]$usr);
         parMfg(mfg=c(1,1), usePar=parList[["postTopHist"]]);
         ## For some reason, the title gets cropped unless I run it with outer=TRUE first,
         ## thereafter it isn't cropped, although no par() values change at all! Bug.
         title(main="", outer=TRUE, line=0);
         #printDebug(c("par('mar'): ", cPaste(format(digits=2, par('mar')))));
      }
      par("xpd"=TRUE);
      if (!is.null(subMain)) {
         title(main=main, line=lineMain, font.main=font.main, cex.main=cexMain);
         title(main=subMain, line=lineMain-1, cex.main=cexSub, font.main=font.main);
      } else {
         title(main=main, line=lineMain, font.main=font.main, cex.main=cexMain);
      }
      par("xpd"=FALSE);
   }
   
   
   if (doTopHist) {
      par("mar"=origPar$mar);
      par("plt"=origPar$plt);
      par("usr"=origPar$usr);
   }
   if (length(pointSubset) > 0) {
      #printDebug(c("pointSubset: ", cPaste(pointSubset)));
      #vp1 <- list(pointSubset=nameVector(pointSubset, names(xValues1[pointSubset])),
      #            pointSubsetX=nameVector(xValues1[pointSubset], names(xValues1[pointSubset])),
      #            pointSubsetY=nameVector(yValues1[pointSubset], names(xValues1[pointSubset])),
      #            pointSubsetCex=nameVector(pointSubsetCex, names(xValues1[pointSubset])),
      #            pointSubsetCol=nameVector(pvColors1[pointSubset], names(xValues1[pointSubset])));
      #vp1DF <- data.frame(lapply(vp1, function(i){i}));
      vp1DF <- data.frame(check.names=FALSE, stringsAsFactors=FALSE,
         "pointSubset"=pointSubset,
         "pointSubsetX"=xValues1[pointSubset],
         "pointSubsetY"=yValues1[pointSubset],
         "pointSubsetCex"=pointSubsetCex,
         "pointSubsetCol"=pvColors1[pointSubset]);
      rownames(vp1DF) <- names(pointSubset);
      invisible(list(x=xValues1, y=yValues1, col=pvColors1,
         pointSubsetDF=vp1DF,
         parList=parList,
         labelCoords=labelCoords));
   } else if (length(multiGenesUp) > 0 || length(multiGenesDown) > 0) {
      invisible(list(x=xValues1, y=yValues1, col=pvColors1,
         multiGenesUp=multiGenesUp, multiGenesDown=multiGenesDown,
         parList=parList,
         labelCoords=labelCoords));
   } else {
      invisible(list(x=xValues1, y=yValues1, col=pvColors1,
         parList=parList,
         labelCoords=labelCoords));
   }
}

#' Draw block arrows in plot margins
#' 
#' Draw block arrows in plot margins
#' 
#' This function draws block arrows in plot margins,
#' intended to be used to describe plot axis ranges
#' with an optional label. The driving example is to
#' describe the number of statistical hits shown
#' on a volcano plot.
#' 
#' Run `blockArrowMargin(doExample=TRUE)` to see a visual
#' example.
#' 
#' @examples
#' blockArrowMargin(doExample=TRUE)
#' 
#' @export
blockArrowMargin <- function
(axisPosition="rightAxis",
 arrowPosition="top",
 arrowLabel="",
 arrowDirection="updown",
 xleft=NULL,
 xright=NULL,
 ybottom=NULL,
 ytop=NULL,
 col="#660000FF",
 labelFont=1,
 labelCex=1,
 border="#000000FF",
 xpd=TRUE,
 parUsr=par("usr"),
 arrowWidthPercent=6,
 blockWidthPercent=5,
 arrowLengthPercent=blockWidthPercent*0.6,
 bufferPercent=0.5,
 blankFirst=FALSE,
 arrowLabelColor=NULL,#"#FFFFFFFF",
 arrowLabelBorder="#000000FF",
 doExample=FALSE,
 doBlockGradient=TRUE,
 doShadowText=TRUE,
 gradientDarkFactor=1.5,
 gradientSFactor=1.5,
 verbose=FALSE,
 ...)
{
   ## Purpose is to draw a block arrow, like ones you see in PowerPoint, but using
   ## rect() syntax as if drawing a rectangle.
   ##
   ## Currently the function draws block arrows outside the plot, as if to label
   ## an axis.
   ##
   ## Trim some of the width away to give room for the arrow to be drawn.
   ##
   ## Some examples:
   if (doShadowText) {
      text_fun <- jamba::shadowText;
   } else {
      text_fun <- text;
   }
   if (doExample) {
      jamba::nullPlot();
      blockArrowMargin(axisPosition="rightAxis",
         ybottom=1.52, 
         arrowPosition="top", 
         col="#BB0000FF", 
         arrowLabel="Up-regulated",
         arrowWidthPercent=arrowWidthPercent, 
         arrowLengthPercent=arrowLengthPercent,
         blockWidthPercent=blockWidthPercent, 
         labelCex=labelCex, 
         ...);
      blockArrowMargin(axisPosition="rightAxis", 
         ytop=1.48, 
         arrowPosition="bottom", 
         col="#0000BBFF", 
         arrowLabel="Down-regulated",
         arrowWidthPercent=arrowWidthPercent, 
         arrowLengthPercent=arrowLengthPercent,
         blockWidthPercent=blockWidthPercent, 
         labelCex=labelCex, 
         ...);
      blockArrowMargin(axisPosition="topAxis", 
         xright=1.48, 
         arrowPosition="left", 
         col="#0000BBFF", 
         arrowLabel="Down-regulated",
         arrowWidthPercent=arrowWidthPercent,
         arrowLengthPercent=arrowLengthPercent,
         blockWidthPercent=blockWidthPercent, 
         labelCex=labelCex, 
         ...);
      blockArrowMargin(axisPosition="topAxis", 
         xleft=1.52, 
         arrowPosition="right", 
         col="#BB0000FF", 
         arrowLabel="Up-regulated",
         arrowWidthPercent=arrowWidthPercent, 
         arrowLengthPercent=arrowLengthPercent,
         blockWidthPercent=blockWidthPercent, 
         labelCex=labelCex, 
         ...);
      return(NULL);
   }
   if (length(axisPosition) > 1) {
      retVals <- lapply(axisPosition, function(iAxisPosition) {
         blockArrowMargin(axisPosition=iAxisPosition, 
            arrowPosition=arrowPosition, 
            arrowLabel=arrowLabel,
            arrowDirection=arrowDirection, 
            xleft=xleft, 
            xright=xright, 
            ybottom=ybottom, 
            ytop=ytop,
            col=col, 
            labelFont=1, 
            labelCex=1, 
            border=border, 
            xpd=xpd, 
            parUsr=parUsr,
            arrowWidthPercent=arrowWidthPercent, 
            arrowLengthPercent=arrowLengthPercent,
            blockWidthPercent=blockWidthPercent, 
            bufferPercent=bufferPercent,
            blankFirst=blankFirst,
            arrowLabelColor=arrowLabelColor, 
            arrowLabelBorder=arrowLabelBorder,
            doExample=doExample, 
            doBlockGradient=doBlockGradient, 
            gradientDarkFactor=gradientDarkFactor, 
            gradientSFactor=gradientSFactor, 
            ...);
      });
      invisible(retVals);
   }
   
   if (doBlockGradient) {
      col <- jamba::alpha2col(
         alpha=jamba::col2alpha(col),
         jamba::fixYellow(col));
      if (length(col) == 1) {
         col2 <- jamba::fixYellow(
            jamba::makeColorDarker(
               col,
               darkFactor=gradientDarkFactor, 
               sFactor=gradientSFactor));
         col2 <- jamba::alpha2col(col2,
            alpha=jamba::col2alpha(col));
         if (jamba::igrepHas("topbottom|bottomtop|leftright|rightleft", arrowPosition)) {
            col <- c(col2, col, col, col2);
         } else {
            col <- c(col, col2);
         }
         colGradient <- colorRampPalette(col,
            alpha=TRUE)(35);
      } else {
         colGradient <- colorRampPalette(col,
            alpha=TRUE)(35);
         col2 <- tail(colGradient, 1);
         col <- head(colGradient, 1);
      }
      col1 <- col[1];
      if (length(arrowLabelColor) == 0) {
         arrowLabelColor <- jamba::setTextContrastColor(head(col2, 1));
      }
   }
   if (length(arrowLabelColor) == 0) {
      arrowLabelColor <- jamba::setTextContrastColor(head(col, 1));
   }

   arrowSets <- list("empty"=list(x=NULL, y=NULL));
   plotWidth <- diff(parUsr[1:2]) / diff(par("plt")[1:2]);
   plotHeight <- diff(parUsr[3:4]) / diff(par("plt")[3:4]);
   if (jamba::igrepHas("rightAxis|leftAxis", axisPosition)) {
      if (is.null(ybottom)) {
         ybottom <- parUsr[3];
      }
      if (is.null(ytop)) {
         ytop <- parUsr[4];
      }
      if (jamba::igrepHas("rightAxis", axisPosition)) {
         if (is.null(xleft)) {
            xleft <- parUsr[2] + 
               plotWidth * bufferPercent / 100;
         }
         if (is.null(xright)) {
            xright <- parUsr[2] + 
               plotWidth * blockWidthPercent/100 +
               plotWidth * bufferPercent / 100;
         }
      } else {
         if (is.null(xright)) {
            xright <- parUsr[1] -
               plotWidth * bufferPercent / 100;
         }
         if (is.null(xleft)) {
            xleft <- parUsr[1] -
               plotWidth * blockWidthPercent/100 -
               plotWidth * bufferPercent / 100;
         }
      }
      arrowDirection <- "updown";
      srtLabel <- 90;
      gradientXY <- "y";
   } else if (jamba::igrepHas("topAxis|bottomAxis", axisPosition)) {
      if (jamba::igrepHas("topAxis", axisPosition)) {
         if (is.null(ybottom)) {
            ybottom <- parUsr[4] +
               plotHeight * bufferPercent / 100;
         }
         if (is.null(ytop)) {
            ytop <- parUsr[4] +
               plotHeight * blockWidthPercent/100 +
               plotHeight * bufferPercent / 100;
         }
      } else {
         if (is.null(ytop)) {
            ytop <- parUsr[3] -
               plotHeight * bufferPercent / 100;
         }
         if (is.null(ybottom)) {
            ybottom <- parUsr[3] -
               plotHeight * blockWidthPercent/100 -
               plotHeight * bufferPercent / 100;
         }
      }
      if (is.null(xleft)) {
         xleft <- parUsr[1];
      }
      if (is.null(xright)) {
         xright <- parUsr[2];
      }
      arrowDirection <- "leftright";
      srtLabel <- 0;
      gradientXY <- "x";
   }
   if (blankFirst) {
      #printDebug("blanking the area first");
      polygon(x=c(xleft, xright, xright, xleft),
         y=c(ybottom, ybottom, ytop, ytop), 
         col="white", 
         border="white", 
         xpd=TRUE);
   }
   if (length(jamba::igrep("up|down", arrowDirection)) > 0) {
      ## Trim away the y-coordinates
      arrowWidth <- abs(xright - xleft);
      arrowWidthDiff <- arrowWidth * (arrowWidthPercent / 50);
      if (xright > xleft) {
         xright1 <- xright - arrowWidthDiff;
         xleft1 <- xleft + arrowWidthDiff;
      } else {
         xright1 <- xright + arrowWidthDiff;
         xleft1 <- xleft - arrowWidthDiff;
      }
      arrowLength <- abs(ytop - ybottom);
      arrowLengthDiff <- plotHeight * arrowLengthPercent / 100;
      if (arrowLengthDiff >= arrowLength) {
         arrowLengthDiff <- arrowLength / 3;
      }
      if (ytop > ybottom) {
         if (length(jamba::igrep("top", arrowPosition)) > 0) {
            ytop1 <- ytop - arrowLengthDiff;
            yArrowPoints <- c(ytop1, ytop1, ytop, ytop1, ytop1);
            xArrowPoints <- c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop1, ybottom));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ytop, ytop);
            xArrowPoints <- c(xleft1, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
         if (length(jamba::igrep("bottom", arrowPosition)) > 0) {
            ybottom1 <- ybottom + arrowLengthDiff;
            yArrowPoints <- c(ybottom1, ybottom1, ybottom, ybottom1, ybottom1);
            xArrowPoints <- rev(c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1));
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom1));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ybottom, ybottom);
            xArrowPoints <- c(xright1, xleft1);
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
      } else {
         if (length(jamba::igrep("bottom", arrowPosition)) > 0) {
            ytop1 <- ytop + arrowLengthDiff;
            yArrowPoints <- c(ytop1, ytop1, ytop, ytop1, ytop1);
            xArrowPoints <- rev(c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1));
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop1, ybottom));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ytop, ytop);
            xArrowPoints <- c(xright1, xleft1);
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
         if (length(jamba::igrep("top", arrowPosition)) > 0) {
            ybottom1 <- ybottom - arrowLengthDiff;
            yArrowPoints <- c(ybottom1, ybottom1, ybottom, ybottom1, ybottom1);
            xArrowPoints <- c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom1));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ybottom, ybottom);
            xArrowPoints <- c(xleft1, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
      }
   } else {
      ## Arrow position is left-to-right
      ## Trim away the x-coordinates
      arrowWidth <- abs(ytop - ybottom);
      arrowWidthDiff <- arrowWidth * (arrowWidthPercent / 50);
      if (ytop > ybottom) {
         ytop1 <- ytop - arrowWidthDiff;
         ybottom1 <- ybottom + arrowWidthDiff;
      } else {
         ytop1 <- ytop + arrowWidthDiff;
         ybottom1 <- ybottom - arrowWidthDiff;
      }
      arrowLength <- abs(xright - xleft);
      arrowLengthDiff <- plotWidth * arrowLengthPercent / 100;
      if (arrowLengthDiff >= arrowLength) {
         arrowLengthDiff <- arrowLength / 3;
      }

      if (xright > xleft) {
         if (length(jamba::igrep("right", arrowPosition)) > 0) {
            xright1 <- xright - arrowLengthDiff;
            xArrowPoints <- c(xright1, xright1, xright, xright1, xright1);
            yArrowPoints <- c(ybottom1, ybottom, mean(c(ybottom, ytop)), ytop, ytop1);
            arrowSets <- c(arrowSets, list("right"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright1));
         } else {
            xright1 <- xright;
            xArrowPoints <- c(xright, xright);
            yArrowPoints <- c(ybottom1, ytop1);
            arrowSets <- c(arrowSets, list("right"=list(x=xArrowPoints, y=yArrowPoints)));
            #print(arrowSets);
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
         if (length(jamba::igrep("left", arrowPosition)) > 0) {
            xleft1 <- xleft + arrowLengthDiff;
            xArrowPoints <- c(xleft1, xleft1, xleft, xleft1, xleft1);
            yArrowPoints <- rev(c(ybottom1, ybottom, mean(c(ybottom, ytop)), ytop, ytop1));
            arrowSets <- c(arrowSets, list("left"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft1, xright1));
         } else {
            xleft1 <- xleft;
            xArrowPoints <- c(xleft, xleft);
            yArrowPoints <- rev(c(ybottom1, ytop1));
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft1, xright1));
         }
      } else {
      }
      
   }
   yLabel <- mean(c(ytop, ybottom));
   xLabel <- mean(c(xleft, xright));
   
   allX <- unlist(lapply(jamba::unvigrep("^empty$", names(arrowSets)), function(i){
      j <- arrowSets[[i]];
      j["x"];
   }))
   allY <- unlist(lapply(jamba::unvigrep("^empty$", names(arrowSets)), function(i){
      j <- arrowSets[[i]];
      j["y"];
   }));
   if (doBlockGradient) {
      arrowSides <- sapply(jamba::unvigrep("^empty$", names(arrowSets)), function(i){
         length(arrowSets[[i]][["x"]]) > 2;
      });
      arrowSides <- paste(names(arrowSides)[arrowSides], collapse="");
      if (arrowSides %in% c("bottom", "left")) {
         colGradient <- rev(colGradient);
         col21 <- col1;
         col1 <- col2;
         col2 <- col1;
      }
      arrowBox <- lapply(jamba::nameVector(jamba::unvigrep("^empty$", names(arrowSets))), function(i){
         xi <- arrowSets[[i]][["x"]];
         yi <- arrowSets[[i]][["y"]];
         xI <- unique(c(xi[1], xi[length(xi)]));
         yI <- unique(c(yi[1], yi[length(yi)]));
         list(x=xI, y=yI);
      });
      arrowBoxX <- unique(sort(c(arrowBox[[1]][["x"]], arrowBox[[2]][["x"]])));
      arrowBoxY <- unique(sort(c(arrowBox[[1]][["y"]], arrowBox[[2]][["y"]])));
      
      xpdPar <- par("xpd");
      par("xpd"=TRUE);
      polygon(x=allX, 
         y=allY, 
         col=tail(colGradient,1), 
         border=border, 
         xpd=TRUE);
      gr1 <- gradient_rect(col=colGradient, 
         gradient=gradientXY, 
         xleft=arrowBoxX[1], 
         xright=arrowBoxX[2],
         ybottom=arrowBoxY[1], 
         ytop=arrowBoxY[2], 
         border="#00000000");
      par("xpd"=xpdPar);
      arrowsDrawn1 <- lapply(jamba::nameVector(jamba::vgrep("top|right", names(arrowSets))), function(i){
         as1 <- arrowSets[[i]];
         polygon(x=as1[["x"]], 
            y=as1[["y"]], 
            col=col2, 
            border=NA, 
            xpd=TRUE);
         as1;
      });
      arrowsDrawn2 <- lapply(jamba::nameVector(jamba::vgrep("bottom|left", names(arrowSets))), function(i){
         as1 <- arrowSets[[i]];
         polygon(x=as1[["x"]], 
            y=as1[["y"]], 
            col=col1, 
            border=col1, 
            xpd=TRUE);
         as1;
      });
      polygon(x=allX, 
         y=allY, 
         col=NA, 
         border=border, 
         xpd=TRUE);
   } else {
      polygon(x=allX, 
         y=allY, 
         col=col, 
         border=border, 
         xpd=TRUE);
      arrowsDrawn1 <- NULL;
      arrowsDrawn2 <- NULL;
   }
   
   ## Optionally label the block arrows
   if (!is.null(arrowLabel) && !arrowLabel %in% c(NA, "")) {
      if (verbose) {
         jamba::printDebug("xLabel: ", 
            round(digits=2, xLabel), 
            ",  yLabel: ", 
            round(digits=2, yLabel));
      }
      text_fun(label=arrowLabel,
         x=xLabel,
         y=yLabel,
         srt=srtLabel,
         xpd=TRUE,
         col=arrowLabelColor,
         font=labelFont,
         cex=labelCex,
         #bg=arrowLabelBorder,
         adj=c(0.5,0.5), 
         ...);
   }
   retVals <- list(arrowSets=arrowSets, 
      arrowsDrawn1=arrowsDrawn1, 
      arrowsDrawn2=arrowsDrawn2,
      allX=allX, 
      allY=allY);
   invisible(retVals);
}

#' Rectangle with color gradient fill
#' 
#' Rectangle with color gradient fill
#' 
#' This function was inspired by the `plotrix::gradient.rect()`
#' function in the plotrix R package. The function is
#' simplified here, and requires a vector of colors in `col`.
#' 
#' @param xleft,ybottom,xright,ytop `numeric` vectors indicating
#'    the position of sides of a rectangle, passed to
#'    `graphics::rect()`. Multiple rectangles may be defined.
#' @param col `character` vector of colors used to fill the rectangles.
#' @param gradient `character` string indicating the direction of
#'    color gradient, with two allowed values: `"x"` and `"y"`.
#' @param border `character` value indicating the color of border
#'    around the rectangle.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' jamba::nullPlot(xlim=c(0,5), ylim=c(0,5), xaxt="s", yaxt="s");
#' gradient_rect(xleft=1, 
#'    ybottom=1, 
#'    xright=2.5, 
#'    ytop=2.5, 
#'    col=jamba::getColorRamp("Reds", n=15))
#' gradient_rect(xleft=2.5, 
#'    ybottom=2.5, 
#'    xright=4, 
#'    ytop=4, 
#'    gradient="y",
#'    col=jamba::getColorRamp("Reds", n=15))
#' 
#' @export
gradient_rect <- function
(xleft,
 ybottom,
 xright,
 ytop,
 col,
 gradient="x",
 border=par("fg"),
 ...)
{
   nslices <- length(col)

   nrect <- max(unlist(lapply(list(xleft, ybottom, xright, ytop),
      length)));
   oldxpd <- par(xpd = NA)
   if (nrect > 1) {
      if (length(xleft) < nrect)
         xleft <- rep(xleft, length.out=nrect)
      if (length(ybottom) < nrect)
         ybottom <- rep(ybottom, length.out=nrect)
      if (length(xright) < nrect)
         xright <- rep(xright, length.out=nrect)
      if (length(ytop) < nrect)
         ytop <- rep(ytop, length.out=nrect)
      for (i in 1:nrect) {
         gradient_rect(xleft[i],
            ybottom[i],
            xright[i],
            ytop[i],
            col=col,
            nslices=nslices,
            gradient=gradient,
            border=border,
            ...)
      }
   } else {
      if (gradient == "x") {
         xinc <- (xright - xleft)/nslices;
         xlefts <- seq(xleft,
            xright - xinc,
            length=nslices);
         xrights <- xlefts + xinc;
         rect(xlefts,
            ybottom,
            xrights,
            ytop,
            col=col,
            lty=0);
         rect(xlefts[1],
            ybottom,
            xrights[nslices],
            ytop,
            border=border);
      } else {
         yinc <- (ytop - ybottom)/nslices;
         ybottoms <- seq(ybottom,
            ytop - yinc,
            length=nslices);
         ytops <- ybottoms + yinc;
         rect(xleft,
            ybottoms,
            xright,
            ytops,
            col=col,
            lty=0);
         rect(xleft,
            ybottoms[1],
            xright,
            ytops[nslices],
            border=border);
      }
   }
   par(oldxpd);
   invisible(col);
}


#' Log-scaled axis including transformed P-values
#' 
#' @export
logAxis <- function
(side,
 at=NULL,
 base=2,
 values=FALSE,
 useFcValues=TRUE,
 las=2,
 makeNegative=FALSE,
 doSignedSignificance=FALSE,
 flipSignedSignificanceSign=FALSE,
 bigMark=",",
 prettyN=5,
 digits=2,
 cex.axis=1,
 font.axis=1,
 ...)
{
   ## Purpose is to draw an axis using log scale
   ## Logic borrowed from 'log10' package, but extended
   ## to allow 2-based (or N-based) log transforms
   if (tolower(side) %in% c("x", "bottom")) {
      side <- 1;
   }
   if (tolower(side) %in% c("y", "left")) {
      side <- 2;
   }
   if (tolower(side) %in% c("above", "top")) {
      side <- 3;
   }
   if (tolower(side) %in% c("right")) {
      side <- 4;
   }
   if (is.null(at)) {
      if (side %in% c(1,3)) {
         #at1 <- rmNA(log(pretty(base^par("usr")[1:2], n=prettyN*10), base=base));
         at <- pretty(c(par("usr")[1:2]), n=prettyN);
      } else {
         #at1 <- rmNA(log(pretty(base^par("usr")[3:4], n=prettyN*10), base=base));
         at <- pretty(c(par("usr")[3:4]), n=prettyN);
      }
      #printDebug(c("at: ", paste(at, collapse=", ")));
   }
   sa1 <- lapply(at, function(i){
      b <- as.numeric(base);
      j <- as.numeric(i);
      #doSignedSignificance <<- doSignedSignificance;
      if (doSignedSignificance %in% c(TRUE, "TRUE")) {
         ## For this operation, force not displaying the value itself
         values <- FALSE;
         ## Grab the sign from the exponent
         jSign <- sign(j);
         ## Make exponents always negative
         j1 <- -(abs(j));
         ## Optionally flip the sign
         if (flipSignedSignificanceSign) {
            b <- jSign * abs(b);
         }
      } else if (makeNegative) {
         j1 <- -j;
      } else {
         j1 <- j;
      }
      if (values) {
         xLabels <- b^j1;
         ## For 2^(-2), instead of reporting 0.5, report -2
         if (useFcValues) {
            xLabels[xLabels < 1] <- (-1 / xLabels[xLabels < 1]);
            xLabels <- signif(xLabels, digits=digits);
         }
         if (!is.null(bigMark) && !bigMark %in% c("")) {
            if (abs(xLabels) >= 1) {
               xLabels <- signif(xLabels, digits=digits);
            }
            xLabels <- format(xLabels, big.mark=bigMark, trim=TRUE);
         }
         axis(side=side, at=j, labels=xLabels, las=2, cex.axis=cex.axis, font.axis=font.axis, ...);
      } else {
         if (i == 0) {
            axis(side=side, at=j, labels=1, las=2, cex.axis=cex.axis*1, font.axis=font.axis, ...);
            xLabels <- 1;
         } else {
            b <- as.character(b);
            j1 <- as.character(j1);
            if (font.axis == 1) {
               axis(side=side, at=j, labels=substitute(b^j1), las=2, cex.axis=cex.axis*1, font.axis=font.axis, ...);
            } else {
               str <- paste0('axis(side, at=j, labels=expression(bold("', b, '"^"', j1, '")), las=2, cex.axis=cex.axis*1, font.axis=2)')
               eval(parse(text=str));
            }
            xLabels <- substitute(b^j1);
         }
      }
      list(j=j, j1=j1, b=b, xLabels=xLabels);
   });
   invisible(list(sa1));
}

