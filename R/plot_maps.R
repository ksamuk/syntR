#' Title
#'
#' @param map_df
#' @param map1_chrom_breaks
#' @param map2_chrom_breaks
#' @param pch_val
#' @param cex_val
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_maps <- function(map_df, map1_chrom_breaks, map2_chrom_breaks, pch_val = 16, cex_val = 0.75,...){

  # check for name consistency and change if necessary

  # allow for 'centroid_x' 'centroid_y' naming for positions
  if("centroid_x" %in% names(map_df)){

    names(map_df)[which(names(map_df) == "centroid_x")] <- "map1_posfull"
    names(map_df)[which(names(map_df) == "centroid_y")] <- "map2_posfull"

  }

  # create map plots
  plot(map_df$map1_posfull, map_df$map2_posfull, pch = pch_val,
       cex = cex_val, xlab = "Map 1", ylab = "Map 2",
        xlim = c(0, max(map1_chrom_breaks)),
        ylim = c(0, max(map2_chrom_breaks)),
       xaxs = "i", yaxs = "i", ...)

  # add chromosome boundary lines
  abline(v = map1_chrom_breaks)
  abline(h = map2_chrom_breaks)

  # add chromosome annotations
  # compute midpoints for each unique chromosome and add text to plot at those points

  midpoints1 <- map_df %>%
    group_by(map1_chr) %>%
    summarize(midpoint = (max(map1_posfull) + min(map1_posfull)) / 2)

  midpoints2 <- map_df %>%
    group_by(map2_chr) %>%
    summarize(midpoint = (max(map2_posfull) + min(map2_posfull)) / 2)

  # add chr1 annotations
  text(x = midpoints1$midpoint, y = max(map2_chrom_breaks) * 1.025,
       labels = midpoints1$map1_chr, xpd = TRUE, cex = 0.75)

  # add chr2 annotations
  text(y = midpoints2$midpoint, x = max(map1_chrom_breaks) * 1.025,
       labels = midpoints2$map2_chr, xpd = TRUE, srt = 270, cex = 0.75)

}
