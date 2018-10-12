#' Determine the appropriate number of markers to include in each initial cluster
#'
#' @param clust
#' @param mark_df
#' @param max_cluster_range
#'
#' @return
#' @export
#'
#' @examples
find_k <- function(clust, mark_df, max_cluster_range) {

  test_k <- function(x){

    mark_df$cluster <- cutree(clust, x)

    cent_df <- mark_df %>%
      group_by(cluster) %>%
      summarise(centroid_x_max = max(sp2), centroid_x_min = min(sp2), centroid_y_max = max(sp1), centroid_y_min = min(sp1))

    return(max(c((cent_df$centroid_x_max - cent_df$centroid_x_min), (cent_df$centroid_y_max - cent_df$centroid_y_min))))

  }

  if (nrow(mark_df) > 1000){

    ks_to_test <- seq(500, nrow(mark_df), by = 500)
    approx_k <- ks_to_test[which.min(abs(unlist(lapply(ks_to_test, test_k)) - max_cluster_range))]
    ks_to_test <- pmin(pmax(seq(approx_k - 500, approx_k + 500, by = 50), 1), nrow(mark_df))
    approx_k <- ks_to_test[which.min(abs(unlist(lapply(ks_to_test, test_k)) - max_cluster_range))]
    ks_to_test <- pmin(pmax(seq(approx_k - 50, approx_k + 50, by = 5), 1), nrow(mark_df))
    approx_k <- ks_to_test[which.min(abs(unlist(lapply(ks_to_test, test_k)) - max_cluster_range))]
    ks_to_test <- pmin(pmax(seq(approx_k - 5, approx_k + 5, by = 1), 1), nrow(mark_df))
    k <- ks_to_test[which.min(abs(unlist(lapply(ks_to_test, test_k)) - max_cluster_range))]

  } else {

    ks_to_test <- seq(50, nrow(mark_df), by = 50)
    approx_k <- ks_to_test[which.min(abs(unlist(lapply(ks_to_test, test_k)) - max_cluster_range))]
    ks_to_test <- pmin(pmax(seq(approx_k - 50, approx_k + 50, by = 5), 1), nrow(mark_df))
    approx_k <- ks_to_test[which.min(abs(unlist(lapply(ks_to_test, test_k)) - max_cluster_range))]
    ks_to_test <- pmin(pmax(seq(approx_k - 5, approx_k + 5, by = 1), 1), nrow(mark_df))
    k <- ks_to_test[which.min(abs(unlist(lapply(ks_to_test, test_k)) - max_cluster_range))]

  }

  return(k)

}
