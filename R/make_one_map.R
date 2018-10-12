#' Collate two genetic maps into a single dataframe
#'
#' @param map A dataframe containing marker positions for two species
#' @param space_size The amount of padding to add between chromosomes (in cM)
#' @param order_chr1 A specific order for the markers in species 1 (optional)
#' @param order_chr2 A specific order for the markers in species 2 (optional)
#' @param max_length1 The maximum length of each chromosome in species 1 (optional)
#' @param max_length2 The maximum length of each chromosome in species 2 (optional)
#' @param flip_chrs A vector of chromosomes to flip in orientation (optional)
#'
#' @return
#' @export
#'
#' @examples
make_one_map <- function(map, space_size = 20, order_chr1 = NULL, order_chr2 = NULL, max_length1 = NULL, max_length2 = NULL, flip_chrs = NULL) {

  # Function to place markers onto a sinlge map
  # map is a dataframe with the columns name1, chr1, pos1, name2, chr2, pos2
  # map1 will be on the x-axis
  # space_size is the amount of space to be added between chromosomes (20 cM default)
  # make sure that max_lengths vectors are in the order that chromosomes will be plotted


  # make sure that map columns are the correct type
  # ? need to include?
  map$chr1 <- as.factor(map$chr1)
  map$chr2 <- as.factor(map$chr2)
  map$pos1 <- as.numeric(map$pos1)
  map$pos2 <- as.numeric(map$pos2)

  # determine the number of chromosomes in each map
  chr_number1 <- length(levels(map$chr1))
  chr_number2 <- length(levels(map$chr2))

  # find max length of each chromosome
  totals1 <- map %>% group_by(chr1) %>% summarize(max = max(pos1))
  totals2 <- map %>% group_by(chr2) %>% summarize(max = max(pos2))

  # if specific order is given change the order of chrosomes
  if (!is.null(order_chr1)) { totals1 <- slice(totals1, match(order_chr1, chr1)) }
  if (!is.null(order_chr2)) { totals2 <- slice(totals2, match(order_chr2, chr2)) }

  # if max length given change chromosome lengths
  if (!is.null(max_length1)) { totals1 <- mutate(totals1, max = max_length1) }
  if (!is.null(max_length2)) { totals2 <- mutate(totals2, max = max_length2) }

  # calculate distance to add to each chromosome
  totals1 <- mutate(totals1, order = 1:chr_number1)
  cumsum1 <- cumsum(totals1$max) + totals1$order * space_size
  totals1$dist2add <- c(0, cumsum1)[1:chr_number1]
  totals2 <- mutate(totals2, order = 1:chr_number2)
  cumsum2 <- cumsum(totals2$max) + totals2$order * space_size
  totals2$dist2add <- c(0, cumsum2)[1:chr_number2]

  # if their are chrs to flip, flip them
  if (!is.null(flip_chrs)) {

    map$pos2flip <- map$pos2

    for(i in flip_chrs) {

      map[map$chr2==i,"pos2flip"] <- as.numeric(totals2[totals2$chr2 == i, "max"]) - map[map$chr2 == i,"pos2"]

    }

    map <- map[,c("name1", "chr1", "pos1", "name2", "chr2", "pos2flip")]
    names(map) <- c("name1", "chr1", "pos1", "name2", "chr2", "pos2")

    }

  # find and replace chr names with distance to be addeded
  find_replace_map1 = setNames(totals1$dist2add, totals1$chr1)
  t1 <- ordered(map$chr1, levels = names(find_replace_map1))
  add1 <- find_replace_map1[t1]
  find.replace.map2 = setNames(totals2$dist2add, totals2$chr2)
  t2 <- ordered(map$chr2, levels = names(find.replace.map2))
  add2 <- find.replace.map2[t2]

  # add extra distance to each postition
  map$pos1full <- map$pos1 + add1
  map$pos2full <- map$pos2 + add2

  # determine the midpoint between each chromosome
  breaks1 <- c(0 - space_size/2, cumsum1 - space_size/2)
  breaks2 <- c(0 - space_size/2, cumsum2 - space_size/2)

  # make a vector of the chr orders
  sp1_chr_order <- totals1 %>% arrange(order) %>% pull(chr1)
  sp2_chr_order <- totals2 %>% arrange(order) %>% pull(chr2)

  full_map_list <- list(map, breaks1, breaks2, sp1_chr_order, sp2_chr_order)

  return(full_map_list)

}







