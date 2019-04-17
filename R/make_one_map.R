#' Consolidate and prepare map data for synteny block detection
#'
#' @param map data frame, contains marker order and chromosome assignments for markers in both maps
#' @param space_size numeric, the amount of space to pad markers between linkage groups
#' @param map1_chr_order numeric, optional: the numeric order in which to display markers in map1
#' @param map2_chr_order numeric, optional: the numeric order in which to display markers in map2
#' @param map1_max_chr_lengths numeric, optional: values of the total lengths of each linkage group in map1. useful when markers do not span known extent of linkage groups.
#' @param map2_max_chr_lengths numeric, optional: values of the total lengths of each linkage group in map2. useful when markers do not span known extent of linkage groups.
#' @param flip_chrs logical, optional: a vector of TRUE/FALSE indicating which linkage groups should be flipped (for display purposes only)
#'
#' @return a 'map list': a list of length 5 containing consolidated/flipped/reordered markers for both maps
#' @export
#'
#' @examples
make_one_map <- function(map, space_size = 20, map1_chr_order = NULL, map2_chr_order = NULL, map1_max_chr_lengths = NULL, map2_max_chr_lengths = NULL, flip_chrs = NULL) {

  # make sure that map columns are the correct type
  map$map1_chr <- as.factor(map$map1_chr)
  map$map2_chr <- as.factor(map$map2_chr)
  map$map1_pos <- as.numeric(map$map1_pos)
  map$map2_pos <- as.numeric(map$map2_pos)

  # determine the number of chromosomes in each map
  map1_chr_number <- length(levels(map$map1_chr))
  map2_chr_number <- length(levels(map$map2_chr))

  # determine the max length of each chromosome
  map1_totals <- map %>% group_by(map1_chr) %>% summarize(map1_max_chr_lengths = max(map1_pos))
  map2_totals <- map %>% group_by(map2_chr) %>% summarize(map2_max_chr_lengths = max(map2_pos))

  # if a specific order is given reorder the chromosomes
  if (!is.null(map1_chr_order)) { map1_totals <- slice(map1_totals, match(map1_chr_order, map1_chr)) }
  if (!is.null(map2_chr_order)) { map2_totals <- slice(map2_totals, match(map2_chr_order, map2_chr)) }

  # if max chromosome lengths were provided, overwrite the automatically calculated values
  if (!is.null(map1_max_chr_lengths)) {

    map1_totals <- map1_totals %>%
      select(map1_chr) %>%
      left_join(map1_max_chr_lengths, by = "map1_chr")

    }
  if (!is.null(map2_max_chr_lengths )) {

    map2_totals <- map2_totals %>%
      select(map2_chr) %>%
      left_join(map2_max_chr_lengths, by = "map2_chr")

    }

  # calculate padding and cumulative chromosome lengths
  # to add to each chromosome for (later) repositioning
  map1_totals <- mutate(map1_totals, order = 1:map1_chr_number)
  map1_cumsum <- cumsum(map1_totals$map1_max_chr_lengths) + map1_totals$order * space_size
  map1_totals$map1_dist2add <- c(0, map1_cumsum)[1:map1_chr_number]

  map2_totals <- mutate(map2_totals, order = 1:map2_chr_number)
  map2_cumsum <- cumsum(map2_totals$map2_max_chr_lengths) + map2_totals$order * space_size
  map2_totals$map2_dist2add <- c(0, map2_cumsum)[1:map2_chr_number]

  # if their are chrs to flip, flip them
  # only applied to one set (map2)

  if (!is.null(flip_chrs)) {

    map$pos2flip <- map$map2_pos

    for(i in flip_chrs) {

      map[map$map2_chr==i,"pos2flip"] <- as.numeric(map2_totals[map2_totals$map2_chr == i, "map2_max_chr_lengths"]) - map[map$map2_chr == i,"map2_pos"]

    }

    map <- map[,c("map1_name", "map1_chr", "map1_pos", "map2_name", "map2_chr", "pos2flip")]
    names(map) <- c("map1_name", "map1_chr", "map1_pos", "map2_name", "map2_chr", "map2_pos")

    }

  # add extra distance to each position in map
  map <- left_join(map, map1_totals, by = "map1_chr") %>%
    left_join(., map2_totals, by = "map2_chr") %>%
    mutate(map1_posfull = map1_pos + map1_dist2add, map2_posfull = map2_pos + map2_dist2add) %>%
    select(map1_name, map1_chr, map1_pos, map1_posfull, map2_name, map2_chr, map2_pos, map2_posfull)

  # determine the midpoint between each chromosome
  map1_breaks <- c(0 - space_size/2, map1_cumsum - space_size/2)
  map2_breaks <- c(0 - space_size/2, map2_cumsum - space_size/2)

  # make a vector of the chr orders
  map1_chr_order <- map1_totals %>% arrange(order) %>% pull(map1_chr)
  map2_chr_order <- map2_totals %>% arrange(order) %>% pull(map2_chr)

  full_map_list <- list(map, map1_breaks, map2_breaks, map1_chr_order, map2_chr_order)

  return(full_map_list)

}







