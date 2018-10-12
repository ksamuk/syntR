#' Estimate mapping error rate from empirical data
#'
#' @param map_list
#' @param linear_regions
#'
#' @return
#' @export
#'
#' @examples
estimate_error_rate_1 <- function(map_list, linear_regions){

  # compute residuals of a linear fit for *each* chromosome
  # uses rlm for robustness to outliers

  chr_resids <- list()

  for (i in 1:nrow(linear_regions)){

    # filter for target region
    chr_df <- map_list[[1]] %>% filter(pos1full>linear_regions[i,1] & pos1full<linear_regions[i,2]) %>%
      filter(pos2full>linear_regions[i,3] & pos2full<linear_regions[i,4])

    # fit a robust linear model (warnings are suppressed)
    suppressWarnings(chr_lm <- rlm(pos2full ~ pos1full, data = chr_df))

    # plot fits
    plot(chr_df$pos2full ~ chr_df$pos1full)
    abline(chr_lm, col = "red")

    # appends the residuals to the residuals list
    chr_resids[[i]] <- residuals(chr_lm)

  }

  # combine the individual residual vectors
  total_resids <- unlist(chr_resids)

  # output the median absolute deviation (a robust metric of SD)
  # for the combined residuals
  sd_est <- mad(total_resids)

  # adjust sd estimate
  return(list(0.002334 + 0.699239*sd_est, sd_est))

}
