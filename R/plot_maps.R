#' Plot a two-dimension representation of two genetic maps
#'
#' @param map_df
#' @param sp1_chrom_breaks
#' @param sp2_chrom_breaks
#' @param pch_val
#' @param cex_val
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_maps <- function(map_df, sp1_chrom_breaks, sp2_chrom_breaks, pch_val = 16, cex_val = 0.75,...){

  # allow for 'sp1' 'sp2' naming for positions
  if("sp1" %in% names(map_df)){

    names(map_df)[which(names(map_df) == "sp1")] <- "pos1full"
    names(map_df)[which(names(map_df) == "sp2")] <- "pos2full"

  }

  # allow for 'centroid_x' 'centroid_y' naming for positions
  if("centroid_x" %in% names(map_df)){

    names(map_df)[which(names(map_df) == "centroid_x")] <- "pos1full"
    names(map_df)[which(names(map_df) == "centroid_y")] <- "pos2full"

  }

  plot(map_df$pos1full, map_df$pos2full, pch = pch_val,
       cex = cex_val, xlab = "Species 1", ylab = "Species 2",
        xlim = c(0, max(sp1_chrom_breaks)),
        ylim = c(0, max(sp2_chrom_breaks)),
       xaxs = "i", yaxs = "i", ...)


  # add chromosome boundary lines
  abline(v = sp1_chrom_breaks)
  abline(h = sp2_chrom_breaks)

  # add chromosome annotations
  # compute midpoints for each unique chromosome and add text to plot at those points

  midpoints1 <- map_df %>%
    group_by(chr1) %>%
    summarize(midpoint = (max(pos1full) + min(pos1full)) / 2)

  midpoints2 <- map_df %>%
    group_by(chr2) %>%
    summarize(midpoint = (max(pos2full) + min(pos2full)) / 2)

  # add chr1 annotations
  text(x = midpoints1$midpoint, y = max(sp2_chrom_breaks) * 1.025,
       labels = midpoints1$chr1, xpd = TRUE, cex = 0.75)

  # add chr2 annotations
  text(y = midpoints2$midpoint, x = max(sp1_chrom_breaks) * 1.025,
       labels = midpoints2$chr2, xpd = TRUE, srt = 270, cex = 0.75)

}
