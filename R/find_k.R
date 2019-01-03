#' Title
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

  # define a function that outputs the maximum range (in x or y) for the clusters defined by any k

  test_k <- function(k){

    mark_df$cluster <- cutree(clust, k)

    cent_df <- mark_df %>%
      group_by(cluster) %>%
      summarise(centroid_x_max = max(map2_posfull), centroid_x_min = min(map2_posfull), centroid_y_max = max(map1_posfull), centroid_y_min = min(map1_posfull))

    return(max(c((cent_df$centroid_x_max - cent_df$centroid_x_min), (cent_df$centroid_y_max - cent_df$centroid_y_min))))

  }

  # find the maximum cluster range for k values that span the number of markers in the maps and
  # identify the smallest k value that yeilds a cluster range that is below the maximum cut-off set by max_cluster_range
  # (small k values means more markers in each cluster)
  # this is done iteratively, starting with k values that are far apart
  # if there are many markers (>1000), the initially tested k values start even farther apart

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
