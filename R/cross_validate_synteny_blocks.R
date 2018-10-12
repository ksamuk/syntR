#' Perform cross validation of synteny blocks via subsampling
#'
#' @param map_list
#' @param max_cluster_range
#' @param max_nn_dist
#' @param min_block_size
#' @param plots
#' @param buffer_cM
#' @param n_sims
#'
#' @return
#' @export
#'
#' @examples
cross_validate_synteny_blocks <- function(map_list, max_cluster_range, max_nn_dist = 15, min_block_size = 3, plots = FALSE, buffer_cM = 1, n_sims = 100) {

  prop_overlap <- vector()
  prop_common_breaks <- vector()

  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

  for (n in 1:n_sims) {

    # find synteny blocks in subsets of the dat

    # subset A
    map_list_subA <- map_list
    sampled_rows <- sample(1:nrow(map_list[[1]]), nrow(map_list[[1]])/2, replace = FALSE)
    map_list_subA[[1]] <- slice(map_list[[1]], sampled_rows)
    map_subA_all_dfs <- find_synteny_blocks(map_list_subA, max_cluster_range = max_cluster_range, max_nn_dist = max_nn_dist, min_block_size = min_block_size, plots = plots)

    # subset B
    map_list_subB <- map_list
    map_list_subB[[1]] <- setdiff(map_list[[1]], map_list_subA[[1]])
    map_subB_all_dfs <- find_synteny_blocks(map_list_subB, max_cluster_range = max_cluster_range, max_nn_dist = max_nn_dist, min_block_size = min_block_size, plots = plots)

    # quanitfy amount of interval overlap between subsets
    overlaping_intervals <- Intervals(rbind(interval_intersection(Intervals(map_subA_all_dfs[[3]]), Intervals(map_subB_all_dfs[[3]])),
                                            interval_intersection(Intervals(map_subA_all_dfs[[4]]), Intervals(map_subB_all_dfs[[4]]))))
    total_intervals <- Intervals(rbind(interval_union(Intervals(map_subA_all_dfs[[3]]), Intervals(map_subB_all_dfs[[3]])),
                                       interval_union(Intervals(map_subA_all_dfs[[4]]), Intervals(map_subB_all_dfs[[4]]))))
    prop_overlap <- c(prop_overlap, sum(size(overlaping_intervals)) / sum(size(total_intervals)))

    # make slightly more permissive breakpoints
    map_subA_all_dfs[[3]][,1] <- map_subA_all_dfs[[3]][,1] - buffer_cM
    map_subA_all_dfs[[3]][,2] <- map_subA_all_dfs[[3]][,2] + buffer_cM
    map_subA_all_dfs[[4]][,1] <- map_subA_all_dfs[[4]][,1] - buffer_cM
    map_subA_all_dfs[[4]][,2] <- map_subA_all_dfs[[4]][,2] + buffer_cM
    map_subB_all_dfs[[3]][,1] <- map_subB_all_dfs[[3]][,1] - buffer_cM
    map_subB_all_dfs[[3]][,2] <- map_subB_all_dfs[[3]][,2] + buffer_cM
    map_subB_all_dfs[[4]][,1] <- map_subB_all_dfs[[4]][,1] - buffer_cM
    map_subB_all_dfs[[4]][,2] <- map_subB_all_dfs[[4]][,2] + buffer_cM

    # count fraction of invervals with some overlap between subsets
    x_overlap <- interval_overlap(Intervals(map_subA_all_dfs[[3]]), Intervals(map_subB_all_dfs[[3]]))
    y_overlap <- interval_overlap(Intervals(map_subB_all_dfs[[4]]), Intervals(map_subA_all_dfs[[4]]))
    overlap_count <- sum(unlist(lapply(x_overlap, length))>0) + sum(unlist(lapply(y_overlap, length))>0)
    prop_common_breaks <- c(prop_common_breaks, overlap_count / (length(x_overlap) + length(y_overlap)))

    # update progess bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, n)

  }

  output <- data.frame(mean_prop_overlap = mean(prop_overlap), mean_prop_common_breaks = mean(prop_common_breaks))

  close(pb)

  return(output)

}
