#' Title
#'
#' @param map_list
#' @param max_cluster_range_list
#' @param max_nn_dist_list
#' @param min_block_size
#'
#' @return
#' @export
#'
#' @examples
test_parameters <- function(map_list, max_cluster_range_list, max_nn_dist_list, min_block_size = 2, progress_bar = FALSE) {

  # run find synteny blocks on each parameter combination and save the summary statistics

  # initialize the data frame and lists
  block_out_summary <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("max_cluster_range", "nn_dist", "num_blocks", "sp1_coverage", "sp2_coverage", "n_outliers"))
  x_breaks_list <- list()
  y_breaks_list <- list()

  # set-up progress bar
  if(progress_bar){

    progress <- 0
    pb <- txtProgressBar(min = 0, max = length(max_cluster_range_list) * length(max_nn_dist_list), style = 3)

  }
  # find synteny blocks for each parameter combination
  for (i in max_cluster_range_list) {

    for (j in max_nn_dist_list) {

      # find synteny blocks
      blocks_list <- find_synteny_blocks(map_list, max_cluster_range = i, max_nn_dist = j, min_block_size)

      # collect data for each run
      y_out <- list(blocks_list[[4]])
      names(y_out) <- paste0("max_cluster_range_", i, "_&_max_nn_dist_", j)
      x_out <- list(blocks_list[[3]])
      names(x_out) <- paste0("max_cluster_range_", i, "_&_max_nn_dist_", j)
      x_breaks_list <- c(x_breaks_list, x_out)
      y_breaks_list <- c(y_breaks_list, y_out)
      block_out_summary <- rbind(block_out_summary, cbind(data.frame(max_cluster_range = i, max_nn_dist = j), blocks_list[[5]]))

      breaks_list <- list(x_breaks_list, y_breaks_list)

      if(progress_bar){

        # update progess bar
        progress <- progress + 1
        setTxtProgressBar(pb, progress)

      }

    }

  }

  if(progress_bar){close(pb)}

  normalize <- function(x){ (x - min(x)) / (max(x) - min(x)) }
  block_out_summary <- mutate(block_out_summary, composite = normalize(map1_coverage) + normalize(map2_coverage) + 1 - normalize(n_outliers))

  return(list(block_out_summary, breaks_list))

}
