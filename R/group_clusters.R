#' Perform cluster re-grouping via the friends-of-friends approach
#'
#' @param neighbour_list
#'
#' @return
#' @export
#'
#' @examples
group_clusters <- function(neighbour_list){

  group_list <- list()

  for(i in 1:length(neighbour_list)){

    # the original group is the marker and its fifth elements
    group_list[[i]] <- c(neighbour_list[[i]]$marker_id, neighbour_list[[i]]$fifth_element)

    # determine the if the group is "complete"
    # i.e. has no further connections to explore
    # (this will only apply to singletons)
    group_complete <- ifelse(length(neighbour_list[[i]]$fifth_element) == 0, TRUE, FALSE)

    # indexes for all the markers in the lists
    # used for convenience below
    marker_indexes <- lapply(neighbour_list, function(x) x$marker_id)

    while(group_complete == FALSE){

      #the neighbour_list indexes of the markers in the current group
      group_indexes <- which(marker_indexes %in% group_list[[i]])

      # extract the unique fifth elements of (all) the markers in the current group
      new_markers <- lapply(group_indexes, function(x) neighbour_list[[x]]$fifth_element) %>% unlist %>% unique

      # if there are any markers in the new list that don't exist in the old one, add them and continue to iterate
      # otherwise end the iteration and sort the final vector
      if(any(!(new_markers %in% group_list[[i]]))){

        group_list[[i]] <- c(group_list[[i]], new_markers) %>% unique

      } else{

        group_complete <- TRUE
        group_list[[i]] <- sort(group_list[[i]])

      }

    }

  }

  return(group_list)

}
