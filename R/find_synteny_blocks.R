#' Indentify synteny blocks via comparison of two genetics maps
#'
#' @param map_list
#' @param max_cluster_range
#' @param max_nn_dist
#' @param min_block_size
#' @param plots
#'
#' @return
#' @export
#'
#' @examples
find_synteny_blocks <- function(map_list, max_cluster_range, max_nn_dist = 15, min_block_size = 2, plots = FALSE) {

  map <- map_list[[1]]
  sp1_chrom_breaks <- map_list[[2]]
  sp2_chrom_breaks <- map_list[[3]]

  # define initial clusters
  mark_df <- data.frame(sp1 = map$pos1full, sp2 = map$pos2full)
  clust <- hclust(dist(mark_df))

  # choose k based on a max range in x or y for a cluster of markers
  k <- find_k(clust, mark_df, max_cluster_range)
  mark_df$cluster <- cutree(clust, k)

  # plot initial clusters (optional)
  if (plots == TRUE){

    # create a (scrambled) list of colors for plotting clusters
    cluster_cols <- sample(viridis(length(unique(mark_df$cluster))))

    # plot the initial cluster assignments
    mark_df %>%
      rename(pos1full = sp1, pos2full = sp2) %>%
      left_join(map, ., by = c("pos1full", "pos2full")) %>%
      plot_maps(sp1_chrom_breaks, sp2_chrom_breaks,
                col = cluster_cols[.$cluster], main = "Initial cluster assignments")

  }

  # determine the centroids for each cluster
  clust_df <- mark_df %>%
    group_by(cluster) %>%
    summarise(centroid_x = mean(sp1), centroid_y = mean(sp2))

  # find nearest neighbour distances for each cluster
  dist_mat <- dist(data.frame(x = clust_df$centroid_x, y = clust_df$centroid_y), method = "maximum")
  dist_df <- melt(as.matrix(dist_mat), varnames = c("cluster", "cluster2"))
  dist_df <- dist_df[dist_df$cluster != dist_df$cluster2,]
  dist_df <- dist_df %>%
    group_by(cluster) %>%
    summarise(distance = min(value))
  clust_df <- clust_df %>%
    left_join(dist_df, by = "cluster")

  # plot clusters that will be dropped (optional)
  if (plots == TRUE){

    # first plot a blank version of mark df
    mark_df %>%
      rename(pos1full = sp1, pos2full = sp2) %>%
      left_join(map, ., by = c("pos1full", "pos2full")) %>%
      plot_maps(sp1_chrom_breaks, sp2_chrom_breaks,
                col = 0, main = "Outlier clusters flagged for removal")

    # add in dropped/retained clusters
    points(clust_df$centroid_x, clust_df$centroid_y,
           cex = 0.75, col = as.numeric(clust_df$distance >= max_nn_dist)+1, pch = 16)

  }

  # filter clusters based on nearest neighbour distance
  clust_df <- clust_df %>% filter(distance <= max_nn_dist)

  # group clusters into blocks based on ranks
  clust_df <- find_initial_blocks(clust_df)

  # sequentially remove blocks made up of too few clusters and group clusters again
  if (min_block_size > 1) {

    for (i in 2:min_block_size) {

      clust_df <- clust_df[!(clust_df$block %in% names(which(table(clust_df$block) < i))),c("cluster", "centroid_x", "centroid_y")]
      clust_df <- find_initial_blocks(clust_df)

    }

  }

  # adjust blocks based on chromsome breaks and overlaps
  adjusted_block_list <- adjust_blocks(map_list, mark_df, clust_df)
  mark_df <- adjusted_block_list[[1]]
  synteny_blocks_df <- adjusted_block_list[[2]]

  # find breakpoints in sp1
  sp1_chrom_breaks <- map_list[[2]]
  synteny_blocks_df <- synteny_blocks_df %>% arrange(x_start)
  sp1_breaks <- Intervals(cbind(head(synteny_blocks_df$x_end, -1), tail(synteny_blocks_df$x_start, -1)))
  sp1_chrom_break_indices <- unlist(interval_overlap(sp1_chrom_breaks, sp1_breaks))
  if (length(sp1_chrom_break_indices) == 0) { sp1_breaks <- as.data.frame(sp1_breaks) } else {
    sp1_breaks <- as.data.frame(sp1_breaks[-sp1_chrom_break_indices,])
  }
  names(sp1_breaks) <- c("start", "end")

  # find breakpoints in sp2
  sp2_chrom_breaks <- map_list[[3]]
  synteny_blocks_df <- synteny_blocks_df %>% arrange(y_start)
  sp2_breaks <- Intervals(cbind(head(synteny_blocks_df$y_end, -1), tail(synteny_blocks_df$y_start, -1)))
  sp2_chrom_break_indices <- unlist(interval_overlap(sp2_chrom_breaks, sp2_breaks))
  if (length(sp2_chrom_break_indices) == 0) { sp2_breaks <- as.data.frame(sp2_breaks) } else {
    sp2_breaks <- as.data.frame(sp2_breaks[-sp2_chrom_break_indices,])
  }
  names(sp2_breaks) <- c("start", "end")

  # add original data and blocks to mark_df
  mark_df <- mark_df %>%
    rename(pos1full = sp1, pos2full = sp2) %>%
    left_join(map, ., by = c("pos1full", "pos2full"))

  # put synteny blocks onto individual chromosomes scale
  synteny_blocks_df$sp1_chr_pos <- sapply(synteny_blocks_df$x_start, function(x) max(which(x >= sp1_chrom_breaks)))
  synteny_blocks_df$chr1 <- sapply(synteny_blocks_df$sp1_chr_pos, function(x) map_list[[4]][x])
  synteny_blocks_df$chr1_added <- sapply(synteny_blocks_df$sp1_chr_pos, function(x) sp1_chrom_breaks[x] + abs(sp1_chrom_breaks[1]))
  synteny_blocks_df$sp2_chr_pos <- sapply(synteny_blocks_df$y_start, function(x) max(which(x >= sp2_chrom_breaks)))
  synteny_blocks_df$chr2 <- sapply(synteny_blocks_df$sp2_chr_pos, function(x) map_list[[5]][x])
  synteny_blocks_df$chr2_added <- sapply(synteny_blocks_df$sp2_chr_pos, function(x) sp2_chrom_breaks[x] + abs(sp1_chrom_breaks[1]))
  synteny_blocks_df <- synteny_blocks_df %>%
    mutate(chr1_start = x_start - chr1_added, chr1_end = x_end - chr1_added, chr2_start = y_start - chr2_added, chr2_end = y_end - chr2_added) %>%
    select(block, chr1, chr1_start, chr1_end, chr2, chr2_start, chr2_end)

  # classify synteny blocks based on slope

  # determine slopes and pvals
  # uses ranks to avoid problems with non-linear regions
  slope_df <- mark_df %>%
    group_by(block) %>%
    mutate(pos1 = rank(pos1), pos2 = rank(pos2)) %>%
    do(block_lm_slope = suppressWarnings(summary(lm(pos1 ~ pos2, data = ., na.action = "na.exclude")))$coefficients[,1][2],
       block_lm_pval = suppressWarnings(summary(lm(pos1 ~ pos2, data = ., na.action = "na.exclude")))$coefficients[,4][2])  %>%
    data.frame

  # categorize blocks by slope and pval
  # uses the original, rather than adjusted blocks
  # (renamed variously to facilitate joining)
  slope_df <- slope_df %>%
    mutate(orientation = case_when(block_lm_slope < 0 & block_lm_pval <= 0.05 ~ "negative",
                                   block_lm_slope > 0 & block_lm_pval <= 0.05 ~ "positive",
                                   TRUE ~ "unclassified")) %>%
    mutate(orientation = ifelse(is.na(block), "outlier", orientation)) %>%
    select(block, orientation) %>%
    rename(final_block = block) %>%
    mutate(final_block = as.integer(final_block))

  # join classifications into the mark_df and synteny_blocks_df
  mark_df <- left_join(mark_df, slope_df, by = "final_block")
  synteny_blocks_df <- left_join(synteny_blocks_df, slope_df %>% rename(block = final_block) %>% mutate(block = as.factor(block)), by = "block")

  # plot original synteny blocks (optional)
  if (plots == TRUE){

    block_cols <- c("#cccccc", sample(viridis(length(unique(mark_df$block)))))

    mark_df %>%
      mutate(block = ifelse(is.na(block), 0, block)) %>%
      plot_maps(sp1_chrom_breaks, sp2_chrom_breaks,  col = block_cols[.$block + 1],
                main = "Initial synteny block assignments", pch = 16,  cex_val = 0.75)

  }

  # plot final synteny blocks (optional)
  if (plots == TRUE){

    mark_df %>%
      mutate(final_block = ifelse(is.na(final_block), 0, final_block)) %>%
      plot_maps(sp1_chrom_breaks, sp2_chrom_breaks, col = block_cols[.$final_block + 1],
                main = "Final synteny block assignments",  cex_val = 0.75)

  }

  # plot final block orientations
  if (plots == TRUE){

    class_cols <- c("blue", "grey", "red", "black")

    mark_df %>%
      plot_maps(map_list[[2]], map_list[[3]], col = class_cols[as.numeric(as.factor(.$orientation))],
                main = "Synteny block orientation",  cex_val = 0.75)

  }

  # calculate summary stats
  summary_stats <- data.frame(num_blocks = nrow(synteny_blocks_df)) %>%
    mutate(sp1_coverage = sum(synteny_blocks_df$chr1_end - synteny_blocks_df$chr1_start), sp2_coverage = sum(synteny_blocks_df$chr2_end - synteny_blocks_df$chr2_start))
  summary_stats$n_outliers <- sum(is.na(mark_df$final_block))

  # output final list
  data_frame_list <- list(marker_df = mark_df, synteny_blocks_df = synteny_blocks_df,
                          sp1_breaks = sp1_breaks, sp2_breaks = sp2_breaks,
                          summary_stats = summary_stats)

  return(data_frame_list)

}
