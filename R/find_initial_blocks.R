#' Title
#'
#' @param clust_df
#'
#' @return
#' @export
#'
#' @examples
find_initial_blocks <- function(clust_df) {

  # rank clusters
  clust_df$x_rank <- dense_rank(rank(clust_df$centroid_x, ties.method = "max"))
  clust_df$y_rank <- dense_rank(rank(clust_df$centroid_y, ties.method = "max"))

  # build the list of neighbour connections
  neighbour_list <- lapply(clust_df$cluster, build_neighbour_rank_lists, rank_df = clust_df)

  # groups clusters based on their neighbours
  block_list <- group_clusters(neighbour_list) %>% unique

  # add block ids to clust_df
  add_group_id <- function(index){

    x <- block_list[index]
    x$id <- index
    x <- data.frame(block = x$id, cluster = x[[1]])

  }

  # expand the list in a dataframe for matching to data
  group_df <- lapply(1:length(block_list), add_group_id) %>% bind_rows

  # add block ids to centroid dataframe
  clust_df <- clust_df %>%
    left_join(group_df, by = "cluster")

}
