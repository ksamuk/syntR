#' Simulate two genetic maps that differ at a single inversion
#'
#' @param n_positions
#' @param n_markers
#' @param inversion_size
#' @param error_rate_1
#' @param error_rate_2
#' @param bend_factor
#'
#' @return
#' @export
#'
#' @examples
make_sim_map <- function(n_positions, n_markers, inversion_size, error_rate_1, error_rate_2, bend_factor = 1) {

  # error rate 1 represents the standard deviation of a normoral distribution of marker postion errors
  # error rate 2 represents the propotion of markers that are "misplaced"

  # define start and stop
  inversion_start <- 30
  inversion_end <- inversion_start + inversion_size

  # define relative positions in both maps
  map <- data.frame(
    name1 = paste("Pos", formatC(seq(0.01, 100, by = 0.01), width = 6, format = "f", flag = "0", digits = 2), sep = ""),
    chr1 = rep(1, 10000),
    pos1 = seq(0.01, 100, by = 0.01),
    name2 = paste("Pos", formatC(seq(0.01, 100, by = 0.01), width = 6, format = "f", flag = "0", digits = 2), sep = ""),
    chr2 = rep(1, 10000),
    pos2 = c(seq(0.01, inversion_start, by = 0.01), seq(inversion_end, inversion_start + 0.01, by = -0.01), seq(inversion_end + 0.01, 100, by = 0.01))
  )

  # sample map based on the number of positions
  map <- map[sample(1:10000, n_positions, replace = FALSE), ]

  # resample position to get to number of markers
  map <- map[sample(1:n_positions, n_markers, replace = TRUE), ]

  # move markers based on error rate 1
  map$pos1 <- map$pos1 + round(rnorm(n_markers, mean = 0, sd = error_rate_1), digits = 2)
  map$pos2 <- map$pos2 + round(rnorm(n_markers, mean = 0, sd = error_rate_1), digits = 2)

  # move some proportion of markers (based on error rate 2) randomly in the map
  map[sample(n_markers, n_markers*error_rate_2, replace = FALSE),]$pos2 <- round(runif(n_markers*error_rate_2, 0, 100), digits = 2)

  # change the rare negative numbers to 0
  map$pos1[map$pos1 < 0] <- 0
  map$pos2[map$pos2 < 0] <- 0

  # apply a non-linear transformation (1 = no bend, < 1 = +ve bend, > 1 = -ve bend)

  map$pos2 = (map$pos2)^bend_factor

  simulated_map_list <- list(map, c(inversion_start, inversion_end))

  return(simulated_map_list)

}
