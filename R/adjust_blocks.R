#' Adjust synteny blocks to prevent spanning across chromosomes and overlapping
#'
#' @param map_list A map list object
#' @param mark_df A marker dataframe
#' @param clust_df A cluster dataframe
#'
#' @return An adjusted map list object
#' @export
#'
#' @examples
adjust_blocks <- function(map_list, mark_df, clust_df){

  # split blocks based on chromosome breaks
  clust_df$block_adj_x <- NA
  clust_df$block_adj_y <- NA
  clust_df <- lapply(unique(clust_df$block), split_blocks, clust_df, map_list)
  clust_df <- bind_rows(clust_df)
  clust_df <- clust_df %>%
    rowwise %>%
    mutate(block = as.character(block)) %>%
    mutate(block = paste(block, block_adj_x, block_adj_y)) %>%
    ungroup %>%
    mutate(block = factor(block) %>% as.numeric %>% as.factor) %>%
    select(-block_adj_x, -block_adj_y)

  # add blocks to the mark_df
  mark_df <- clust_df %>%
    select(cluster, block) %>%
    left_join(mark_df, . , by = "cluster")

  # summarize initial synteny blocks
  synteny_blocks_df <- mark_df %>%
    filter(!is.na(block)) %>%
    group_by(block) %>%
    summarise(x_start = min(sp1), x_end = max(sp1), y_start = min(sp2), y_end = max(sp2))

  # drop completely subsumed blocks

  synteny_blocks_df$x_subsummed <- rep(NA, nrow(synteny_blocks_df))
  synteny_blocks_df$y_subsummed <- rep(NA, nrow(synteny_blocks_df))

  for (i in 1:nrow(synteny_blocks_df)){

    other_synteny_blocks <- synteny_blocks_df[-i,]

    x_subsummed <- rep(NA, nrow(other_synteny_blocks))
    y_subsummed <- rep(NA, nrow(other_synteny_blocks))

    for (j in 1:nrow(other_synteny_blocks)){

      x_subsummed[j] <- (other_synteny_blocks$x_start[j] <= synteny_blocks_df$x_start[i]) &
         (other_synteny_blocks$x_end[j] >= synteny_blocks_df$x_end[i])

      y_subsummed[j] <-(other_synteny_blocks$y_start[j] <= synteny_blocks_df$y_start[i]) &
         (other_synteny_blocks$y_end[j] >= synteny_blocks_df$y_end[i])

    }

    synteny_blocks_df$x_subsummed[i] <- any(x_subsummed, na.rm= TRUE)
    synteny_blocks_df$y_subsummed[i] <- any(y_subsummed, na.rm= TRUE)

  }

  # filter subsummed blocks (if any)
  synteny_blocks_df <- synteny_blocks_df %>%
    filter(!(x_subsummed) & !(y_subsummed))

  # adjust synteny blocks starts and ends so that they don't overlap in x
  synteny_blocks_df <- synteny_blocks_df %>%
    arrange(x_start, x_end) %>%
    mutate(new_x_start = ifelse(x_start > lag(x_end), x_start, lag(x_end)),
           new_x_end = ifelse(x_end < lead(x_start), x_end, lead(x_start)))
  synteny_blocks_df$new_x_start[1] <- synteny_blocks_df$x_start[1]
  synteny_blocks_df$new_x_end[length(synteny_blocks_df$new_x_end)] <- synteny_blocks_df$x_end[length(synteny_blocks_df$x_end)]

  # adjust synteny blocks starts and ends so that they don't overlap in y
  synteny_blocks_df <- synteny_blocks_df %>%
    arrange(y_start, y_end) %>%
    mutate(new_y_start = ifelse(y_start > lag(y_end), y_start, lag(y_end)),
           new_y_end = ifelse(y_end < lead(y_start), y_end, lead(y_start)))
  synteny_blocks_df$new_y_start[1] <- synteny_blocks_df$y_start[1]
  synteny_blocks_df$new_y_end[length(synteny_blocks_df$new_y_end)] <- synteny_blocks_df$y_end[length(synteny_blocks_df$y_end)]

  # remove synteny blocks without a positive range in x or y
  synteny_blocks_df <- synteny_blocks_df %>%
    mutate(new_x_range = new_x_end - new_x_start, new_y_range = new_y_end - new_y_start)
  synteny_blocks_df <- synteny_blocks_df[synteny_blocks_df$new_x_range > 0 & synteny_blocks_df$new_y_range > 0, c(1, 8:11)]
  names(synteny_blocks_df)[2:5] <- c("x_start", "x_end", "y_start", "y_end")
  synteny_blocks_df <- as.data.frame(synteny_blocks_df)

  # remove markers from blocks that fall outside adjusted synteny blocks
  mark_df <- left_join(mark_df, synteny_blocks_df, by = "block") %>%
    mutate(final_block = ifelse((sp1 < x_start | sp1 > x_end | sp2 < y_start | sp2 > y_end), NA, block)) %>%
    select(sp1, sp2, cluster, block, final_block)

  out_list <- list(mark_df, synteny_blocks_df)

  return(out_list)

}
