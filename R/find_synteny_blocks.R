#' Title
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
find_synteny_blocks <- function(map_list, max_cluster_range, max_nn_dist, min_block_size = 2, plots = FALSE) {

  # pull out map and chromosome information
  # testing git
  map <- map_list[[1]]
  map1_chrom_breaks <- map_list[[2]]
  map2_chrom_breaks <- map_list[[3]]

  # define initial clusters
  mark_df <- map %>% select(map1_posfull, map2_posfull)
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
      left_join(map, ., by = c("map1_posfull", "map2_posfull")) %>%
      plot_maps(map1_chrom_breaks, map2_chrom_breaks,
                col = cluster_cols[.$cluster], main = "Initial cluster assignments")

  }

  # determine the centroids for each cluster
  clust_df <- mark_df %>%
    group_by(cluster) %>%
    summarise(centroid_x = mean(map1_posfull), centroid_y = mean(map2_posfull))

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
      left_join(map, ., by = c("map1_posfull", "map2_posfull")) %>%
      plot_maps(map1_chrom_breaks, map2_chrom_breaks,
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

  # find breakpoints in map1
  map1_chrom_breaks <- map_list[[2]]
  synteny_blocks_df <- synteny_blocks_df %>% arrange(x_start)
  map1_breaks <- Intervals(cbind(head(synteny_blocks_df$x_end, -1), tail(synteny_blocks_df$x_start, -1)))

  # remove chromosome breaks from breaks list
  map1_chrom_break_indices <- unlist(interval_overlap(map1_chrom_breaks, map1_breaks))
  if (length(map1_chrom_break_indices) == 0) { map1_breaks <- as.data.frame(map1_breaks) } else {
    map1_breaks <- as.data.frame(map1_breaks[-map1_chrom_break_indices,])
  }
  names(map1_breaks) <- c("start", "end")

  # find breakpoints in map2
  map2_chrom_breaks <- map_list[[3]]
  synteny_blocks_df <- synteny_blocks_df %>% arrange(y_start)
  map2_breaks <- Intervals(cbind(head(synteny_blocks_df$y_end, -1), tail(synteny_blocks_df$y_start, -1)))

  # remove chromosome breaks from breaks list
  map2_chrom_break_indices <- unlist(interval_overlap(map2_chrom_breaks, map2_breaks))
  if (length(map2_chrom_break_indices) == 0) { map2_breaks <- as.data.frame(map2_breaks) } else {
    map2_breaks <- as.data.frame(map2_breaks[-map2_chrom_break_indices,])
  }
  names(map2_breaks) <- c("start", "end")

  # add original data and blocks to mark_df
  mark_df <- mark_df %>%
    left_join(map, ., by = c("map1_posfull", "map2_posfull"))

  # put synteny blocks onto individual chromosomes scale (i.e. reverse make_one_map)
  synteny_blocks_df$map1_chr_index <- sapply(synteny_blocks_df$x_start, function(x) max(which(x >= map1_chrom_breaks)))
  synteny_blocks_df$map1_chr <- sapply(synteny_blocks_df$map1_chr_index, function(x) map_list[[4]][x])
  synteny_blocks_df$map1_added <- sapply(synteny_blocks_df$map1_chr_index, function(x) map1_chrom_breaks[x] + abs(map1_chrom_breaks[1]))

  synteny_blocks_df$map2_chr_index <- sapply(synteny_blocks_df$y_start, function(x) max(which(x >= map2_chrom_breaks)))
  synteny_blocks_df$map2_chr <- sapply(synteny_blocks_df$map2_chr_index, function(x) map_list[[5]][x])
  synteny_blocks_df$map2_added <- sapply(synteny_blocks_df$map2_chr_index, function(x) map2_chrom_breaks[x] + abs(map1_chrom_breaks[1]))

  synteny_blocks_df <- synteny_blocks_df %>%
    mutate(map1_start = x_start - map1_added, map1_end = x_end - map1_added, map2_start = y_start - map2_added, map2_end = y_end - map2_added) %>%
    select(block, map1_chr, map1_start, map1_end, map2_chr, map2_start, map2_end)

  # classify synteny blocks based on slope

  # determine slopes and pvals
  # uses ranks to avoid problems with non-linear regions
  slope_df <- mark_df %>%
    group_by(block) %>%
    mutate(map1_pos = rank(map1_pos), map2_pos = rank(map2_pos)) %>%
    do(block_lm_slope = suppressWarnings(summary(lm(map1_pos ~ map2_pos, data = ., na.action = "na.exclude")))$coefficients[,1][2],
       block_lm_pval = suppressWarnings(summary(lm(map1_pos ~ map2_pos, data = ., na.action = "na.exclude")))$coefficients[,4][2])  %>%
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

  # plot final synteny blocks (optional)
  if (plots == TRUE){


    # the palatte creation function
    get_col <- viridis::viridis_pal(alpha = 1.0, begin = 0.15, end = 1.0, direction = 1, option = "D")

    # plot the final clusters, with colors determined per linakge group
    mark_df %>%
      group_by(map1_chr, map2_chr) %>%
      mutate(n_groups = length(unique(final_block)) -1 ) %>%
      mutate(block_col = factor(final_block, levels = unique(final_block[order(map1_pos)])) %>% as.numeric()) %>%
      mutate(block_col = ifelse(is.na(final_block), "#cccccc", get_col(unique(n_groups))[block_col])) %>%
      ungroup %>%
      plot_maps(map1_chrom_breaks, map2_chrom_breaks, col = .$block_col,
                main = "Final synteny block assignments",  cex_val = 1)



  }


  # calculate summary stats
  summary_stats <- data.frame(num_blocks = nrow(synteny_blocks_df)) %>%
    mutate(map1_coverage = sum(synteny_blocks_df$map1_end - synteny_blocks_df$map1_start),
           map2_coverage = sum(synteny_blocks_df$map2_end - synteny_blocks_df$map2_start))

  summary_stats$n_outliers <- sum(is.na(mark_df$final_block))

  # output final list
  data_frame_list <- list(marker_df = mark_df, synteny_blocks_df = synteny_blocks_df,
                          map1_breaks = map1_breaks, map2_breaks = map2_breaks,
                          summary_stats = summary_stats)

  return(data_frame_list)

}
