#' Title
#'
#' @param summary_stats
#' @param stats_to_plot
#'
#' @return
#' @export
#'
#' @examples
plot_summary_stats <- function(summary_stats, stats_to_plot) {

  for (i in stats_to_plot) {

    stat_to_plot <- i

    dat_mat <- select(summary_stats, max_cluster_range, max_nn_dist, stat_to_plot) %>%
      spread(., max_nn_dist, stat_to_plot) %>%
      select(-max_cluster_range) %>%
      data.matrix()

    image(dat_mat, col = cm.colors(256), xaxt = "n", yaxt = "n", xlab = "max_cluster_range", ylab = "max_nn_dist", main = stat_to_plot)

    max_cluster_range_list <- unique(summary_stats$max_cluster_range)
    max_nn_dist_list <- unique(summary_stats$max_nn_dist)

    axis(side = 1, at = seq(0, 1, 1/(length(max_cluster_range_list)-1)), labels = max_cluster_range_list)
    axis(side = 2, at = seq(0, 1, 1/(length(max_nn_dist_list)-1)), labels = max_nn_dist_list)

    for (j in 1:nrow(dat_mat)) {

      for (k in 1:ncol(dat_mat)) {

        text(seq(0, 1, 1/(nrow(dat_mat)-1))[j], seq(0, 1, 1/(ncol(dat_mat)-1))[k], round(dat_mat[j,k], digits = 2))

      }

    }

  }

}
