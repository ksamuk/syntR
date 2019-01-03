#' Title
#'
#' @param query_cluster
#' @param rank_df
#'
#' @return
#' @export
#'
#' @examples
build_neighbour_rank_lists <- function(query_cluster, rank_df){

  # determine the ranks of the query cluster
  cluster_sub <- rank_df %>% filter(cluster == query_cluster)
  cluster_x_rank <- cluster_sub$x_rank
  cluster_y_rank <- cluster_sub$y_rank

  # find the clusters that are 1 rank away (+/-)
  x_neighbour <- rank_df$cluster[rank_df$x_rank == (cluster_x_rank + 1) | rank_df$x_rank == (cluster_x_rank -1)]
  y_neighbour <- rank_df$cluster[rank_df$y_rank == (cluster_y_rank + 1) | rank_df$y_rank == (cluster_y_rank -1)]

  # find the clusters that are 2 rank away (+/-)
  x_p1_neighbour <- rank_df$cluster[rank_df$x_rank == (cluster_x_rank + 2) | rank_df$x_rank == (cluster_x_rank -2)]
  y_p1_neighbour <- rank_df$cluster[rank_df$y_rank == (cluster_y_rank + 2) | rank_df$y_rank == (cluster_y_rank -2)]

  # grouping (5th!) element
  # Make a 5th element that lists the numbers shared in [x]&[y], [x+1]&[y], [x]&[y+1]
  x_y_shared <- x_neighbour[x_neighbour %in% y_neighbour]
  x_p1_y_shared <- x_p1_neighbour[x_p1_neighbour %in% y_neighbour]
  x_y_p1_shared <- x_neighbour[x_neighbour %in% y_p1_neighbour]
  fifth_element <- c(x_y_shared, x_p1_y_shared, x_y_p1_shared) %>% unique

  # output named list of neighbours
  list(marker_id = query_cluster, x_neighbour = x_neighbour, y_neighbour = y_neighbour, x_p1_neighbour = x_p1_neighbour, y_p1_neighbour = y_p1_neighbour, fifth_element = fifth_element)

}
