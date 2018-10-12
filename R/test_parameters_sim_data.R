#' Test the performance of tuning parameters on simulated data
#'
#' @param max_cluster_range_list
#' @param max_nn_dist
#' @param min_block_size
#' @param plots
#' @param n_sims
#' @param buffer_cM
#' @param n_positions
#' @param n_markers
#' @param inversion_size
#' @param error_rate_1
#' @param error_rate_2
#'
#' @return
#' @export
#'
#' @examples
test_parameters_sim_data <- function(max_cluster_range_list, max_nn_dist = 15, min_block_size = 2, plots = FALSE, n_sims, buffer_cM, n_positions, n_markers, inversion_size, error_rate_1, error_rate_2){

  false_positives_sums <- vector()
  false_negatives_sums <- vector()

  for (i in 1:length(max_cluster_range_list)) {

    print(c("max_cluster_range", i), quote = FALSE)

    max_cluster_range <- max_cluster_range_list[i]

    false_positives <- vector()
    false_negatives <- vector()

    pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

    for (n in 1:n_sims) {

      # simiulate map comparison with specific inversions and error rates
      sim_map_list <- make_sim_map(n_positions, n_markers, inversion_size, error_rate_1, error_rate_2)
      map <- sim_map_list[[1]]
      true_breaks <- sim_map_list[[2]]

      # find synteny blocks + breakpoints
      map_list <- make_one_map(map, 0, NULL, NULL)
      map <- map_list[[1]]
      x_breaks <- find_synteny_blocks(map_list, max_cluster_range, max_nn_dist, min_block_size, plots)[[3]]

      if (nrow(x_breaks)== 0) {

        false_positives <- c(false_positives, 0)
        false_negatives <- c(false_negatives, 2)

      } else {

        # make slightly more permissive breakpoints
        approx_x_breaks <- x_breaks
        approx_x_breaks[,1] <- x_breaks[,1] - buffer_cM
        approx_x_breaks[,2] <- x_breaks[,2] + buffer_cM

        # determine the number of false positive and false negatives
        fit <- interval_overlap(Intervals(approx_x_breaks), true_breaks)
        false_positives <- c(false_positives, sum(lengths(fit) == 0))
        false_negatives <- c(false_negatives, max((length(true_breaks) - sum(lengths(fit) == 1)), 0))

      }

      # update progess bar
      setTxtProgressBar(pb, n)

    }

    false_positives_sums <- c(false_positives_sums, sum(false_positives))
    false_negatives_sums <- c(false_negatives_sums, sum(false_negatives))

    close(pb)

  }

  mat_list <- list(false_positive_rate = false_positives_sums/n_sims, false_negative_rate = false_negatives_sums/n_sims)

  return(mat_list)

}
