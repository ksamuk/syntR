#' Split and reorder synteny blocks
#'
#' @param block_target character, an id for a synteny block to scan for breaks
#' @param cent_df data frame, the centroid clustered markers in both maps
#' @param map_list list, the map_list object created by find_synteny_blocks
#'
#' @return a dataframe (sub_df) containing adjusted synteny block assignements for the target block
#' @export
#'
#' @examples
split_blocks <- function(block_target, cent_df, map_list){

  # prevent blocks from occuring across chromosomes
  # scan each block and check if a break (map_list[[2]] and map_list[[3]]) exists within it
  # if a break exists, divide the block into two new blocks

  sub_df <- cent_df %>%
    filter(block == block_target)

  # x positions
  # are any breaks in that range of the block?
  x_breaks <- map_list[[2]][which(map_list[[2]] >= min(sub_df$centroid_x) & map_list[[2]] <= max(sub_df$centroid_x))]
  x_breaks <- sort(x_breaks)

  # if there are any breaks in the range, create a new adjusted group based on each break
  # using the magic of findInterval

  if(length(x_breaks >= 1)){

    sub_df$block_adj_x <- findInterval(sub_df$centroid_x, x_breaks)

  }

  # repeat for y positions
  # are any breaks in the range of the block?
  y_breaks <- map_list[[3]][which(map_list[[3]] >= min(sub_df$centroid_y) & map_list[[3]] <= max(sub_df$centroid_y))]
  y_breaks <- sort(y_breaks)

  if(length(y_breaks >= 1)){

    sub_df$block_adj_y <- findInterval(sub_df$centroid_y, y_breaks)

  }

  return(sub_df)

}
