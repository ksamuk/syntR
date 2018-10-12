#' Estimate mapping noise rate from empirical data
#'
#' @param map_list
#' @param colinear_chrs
#'
#' @return
#' @export
#'
#' @examples
estimate_error_rate_2 <- function(map_list, colinear_chrs){

  map <- map_list[[1]]
  total_markers <- map %>% filter(chr1 %in% colinear_chrs$chr1) %>% nrow

  # find max cM for each chromosome in sp2
  sp2_chr_lengths <- map %>% group_by(chr2) %>%
    summarise(max = max(pos2))

  # find sp2 total map length
  sp2_total_length <- sum(sp2_chr_lengths$max)

  # for each colinear pair count outliers (i.e. markers on other chromosomes) &
  # calculate the adjustment needed to account for the outliers missed on the colinear chromosome
  # this is proportional to sp2 chromosome length relative to sp2 total length
  outliers <- data.frame(num = vector(), adjust = vector())
  for (i in 1:nrow(colinear_chrs)){

    outliers[i,1] <- map %>% filter(chr1==as.character(colinear_chrs[i,1])) %>%
      filter(chr2!=as.character(colinear_chrs[i,2])) %>% count
    outliers[i,2] <- sp2_total_length / (sp2_total_length - sp2_chr_lengths[sp2_chr_lengths$chr2==as.character(colinear_chrs[i,2]),]$max)

  }

  total_outs <- sum(outliers$num * outliers$adjust)

  # calculate the number of outliers in a 100 cM window
  outs_100cM <- total_outs * 100 / sp2_total_length

  # calculate the proportion of markers that would be seen as outliers
  # need to account for the markers that fall outside the 100 cM window
  return(outs_100cM / (total_markers - (total_outs - outs_100cM)))

}

