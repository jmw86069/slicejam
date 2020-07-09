
## slicejam_blendcolors.R

#' Blend multiple colors
#' 
#' Blend multiple colors
#' 
#' This function blends multiple colors, including several useful
#' features:
#' 
#' * color wheel red-yellow-blue which supports subtractive color mixing,
#' as seen when mixing paints
#' * accounts for transparency of individual colors
#' 
#' The basic design guide was to meet these expectations:
#' 
#' * red + yellow = orange
#' * blue + yellow = green
#' * red + blue = purple
#' * blue + red + yellow = some brown/gray substance
#' 
#' The default additive color mixing, with red-green-blue colors
#' used in electronic monitors, does not meet these criteria.
#' (In no logical paint mixing exercise would someone expect that
#' mixing red and green would make yellow.)
#' 
#' @param x `character` vector of R colors in hex format.
#' @param preset `character` value indicating the color wheel preset,
#'    passed to `colorjam::h2hwOptions()`.
#' @param lens `numeric` value used to influence the color saturation
#'    after averaging color wheel angles.
#' @param do_plot `logical` indicating whether to depict the color
#'    blend operation using `jamba::showColors()`.
#' @param c_weight `numeric` value used to weight the average color
#'    chroma (saturation) using the mean chroma values of the input
#'    colors. When `c_weight=0` the chroma uses the radius returned
#'    by the mean color wheel angle.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' blend_colors(c("red", "yellow"), do_plot=TRUE)
#' 
#' blend_colors(c("blue", "yellow"), do_plot=TRUE)
#' 
#' blend_colors(c("blue", "red"), do_plot=TRUE)
#' 
#' blend_colors(c("green4", "red"), do_plot=TRUE)
#' 
#' blend_colors(c("green", "dodgerblue"), do_plot=TRUE)
#' 
#' blend_colors(c("red", "yellow", "blue"), do_plot=TRUE)
#' 
#' @export
blend_colors <- function
(x,
 preset=c("ryb", "none", "dichromat", "rgb", "ryb2"),
 lens=0,
 h1=NULL,
 h2=NULL,
 do_plot=FALSE,
 c_weight=0.2,
 ...)
{
   ## 1. Convert colors to ryb
   ## 2. Implement color subtraction
   ##     new_r <- 255 - sqrt( (255 - r_1)^2 + (255 - r_2)^2 )
   ## 3. convert to rgb
   #x <- jamba::nameVector(c("red", "yellow", "blue"));
   preset <- match.arg(preset);
   
   ## handle list input
   if (is.list(x)) {
      x_unique <- unique(x);
      x_match <- match(x, x_unique);
      x_blends <- sapply(x_unique, function(x1){
         blend_colors(x=x1,
            preset=preset,
            lens=lens,
            h1=h1,
            h2=h2,
            do_plot=FALSE,
            c_weight=c_weight);
      });
      x_blend <- x_blends[x_match];
      names(x_blend) <- names(x);
      return(x_blend);
   }

   ## weights are defined by transparency
   x_w <- jamba::col2alpha(x);
   
   x_rgb <- farver::decode_colour(x, to="rgb");
   #jamba::printDebug("x_rgb:");print(x_rgb);
   x_hcl <- farver::decode_colour(x, to="hcl");
   #jamba::printDebug("x_hcl:");print(x_hcl);
   
   ## adjust hue using color wheel preset
   if (length(h1) == 0 || length(h2) == 0) {
      h1h2 <- colorjam::h2hwOptions(preset=preset,
         setOptions="FALSE");
   } else {
      h1h2 <- list(h1=h1, h2=h2);
   }
   h_rgb <- x_hcl[,"h"];
   h_ryb <- colorjam::h2hw(h=h_rgb,
      h1=h1h2$h1,
      h2=h1h2$h2);
   
   ## mean hue angle
   #printDebug("h_ryb:", round(h_ryb));
   h_ryb_mean_v <- mean_angle(h_ryb,
      lens=lens,
      w=x_w);
   h_ryb_mean <- h_ryb_mean_v["deg"];
   #printDebug("h_ryb_mean:", round(h_ryb_mean));

   h_rgb_mean <- colorjam::hw2h(h=h_ryb_mean,
      h1=h1h2$h1,
      h2=h1h2$h2);
   #printDebug("h_rgb_mean:", round(h_rgb_mean));
   
   mean_radius <- weighted.mean(c(1, h_ryb_mean_v["radius2"]),
      w=c(c_weight, 1));
   #printDebug("radius:", round(h_ryb_mean_v["radius"]*10)/10);
   #printDebug("radius2:", round(h_ryb_mean_v["radius2"]*10)/10);
   #printDebug("mean_radius:", round(mean_radius*10)/10);
   x_mean_c <- weighted.mean(x_hcl[,"c"], w=x_w) * mean_radius;
   x_mean_l <- weighted.mean(x_hcl[,"l"], w=x_w);
   new_hcl <- t(as.matrix(c(h=h_rgb_mean,
      c=unname(x_mean_c),
      l=x_mean_l)));
   #printDebug("new_hcl:");print(new_hcl);
   new_col <- farver::encode_colour(new_hcl, from="hcl");
   if (do_plot) {
      jamba::showColors(list(x=x,
         blended=rep(new_col, length(x))));
   }
   return(new_col);

}

#' Calculate the mean angle
#' 
#' Calculate the mean angle
#' 
#' This function takes a vector of angles in degrees (0 to 360 degrees)
#' and returns the mean angle based upon the average of unit vectors.
#' 
#' The function also optionally accomodates weighted mean values,
#' if a vector of weights is supplied as `w`.
#' 
#' Part of the intent of this function is to be used for color blending
#' methods, for example taking the average color hue from a vector of
#' colors. For this purpose, some colors may have varying color saturation
#' and transparency, which are mapped here as weight `w`. Colors which are
#' fully transparent should therefore have weight `w=0` so they do not
#' contribute to the resulting average color hue. Also during color blending
#' operations, the resulting color saturation is adjusted using the `lens`
#' argument, the default `lens=-5` has a tendency to increase intermediate
#' color saturation.
#' 
#' @return `numeric` vector that contains
#'    * `degree` the mean angle in degrees
#'    * `radius` the actual radius based upon mean unit vectors
#'    * `radius2` the adjusted radius using `jamba::warpAroundZero()`
#' 
#' @param x `numeric` vector of angles in degrees
#' @param w `numeric` vector representing weights
#' @param do_plot `logical` indicating whether to create a visual summary plot
#' @param lens `numeric` value passed to `jamba::warpAroundZero()` to adjust
#'    the radius
#' @param ... additional arguments are ignored
#' 
#' @export
mean_angle <- function
(x,
 w=NULL,
 do_plot=FALSE,
 lens=-5,
 ...)
{
   xy <- data.frame(x=sin(jamba::deg2rad(x)),
      y=cos(jamba::deg2rad(x)));
   if (length(w) == 0) {
      w <- 1;
   }
   w <- rep(w,
      length.out=length(x));

   xy_mean <- matrixStats::colWeightedMeans(
      x=as.matrix(xy),
      w=w);
   xy_m <- matrix(ncol=2, byrow=TRUE,
      c(0, 0, xy_mean));
   x_radius <- dist(xy_m);
   
   x_radius2 <- jamba::warpAroundZero(x_radius,
      xCeiling=1,
      lens=lens);
   
   x_deg <- rad2deg(atan2(x=xy_mean["y"], y=xy_mean["x"])) %% 360;
   
   if (do_plot) {
      jamba::nullPlot(xlim=c(-1,1),
         ylim=c(-1,1),
         asp=1,
         doBoxes=FALSE);
      aseq <- seq(from=0, to=360, by=2);
      lines(x=sin(jamba::deg2rad(aseq)),
         y=cos(jamba::deg2rad(aseq)),
         type="l",
         lty="dotted");
      arrows(x0=0,
         y0=0,
         x1=xy$x * w,
         y1=xy$y * w,
         lwd=2,
         angle=30);
      arrows(x0=0,
         y0=0,
         x1=xy$x,
         y1=xy$y,
         lty="dotted",
         lwd=2,
         angle=90);
      arrows(x0=0,
         y0=0,
         x1=sin(jamba::deg2rad(x_deg)) * x_radius,
         y1=cos(jamba::deg2rad(x_deg)) * x_radius,
         lty="solid",
         lwd=2,
         angle=90,
         col="dodgerblue");
      arrows(x0=0,
         y0=0,
         x1=sin(jamba::deg2rad(x_deg)) * x_radius2,
         y1=cos(jamba::deg2rad(x_deg)) * x_radius2,
         lwd=4,
         angle=30,
         col="darkorange1");
   }
   c(deg=unname(x_deg),
      radius=x_radius,
      radius2=x_radius2);
}

