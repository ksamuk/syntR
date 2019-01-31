## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 6,
  fig.width = 6
)



## ----input_example, echo=FALSE, message=FALSE----------------------------

library("syntR")

data(ann_pet_map)
knitr::kable(head(ann_pet_map), caption = "The example input file")


## ----read_map------------------------------------------------------------
# install the package if needed!
# devtools::install_github("ksamuk/syntR")

# load the syntR package
library("syntR")

# load the example marker data
data(ann_pet_map)
as_tibble(ann_pet_map)

# read in list of chromosome lengths for map1
# this is an optional file that defines the 
# maximum length of each chromosome (in cM) if this is known
data(ann_chr_lengths)
ann_chr_lengths

# make_one_map places all the markers on a single scale
# and adds padding between the chromosomes 
# to aid the algorithm and for visualization purposes
map_list <- make_one_map(ann_pet_map, map1_max_chr_lengths = ann_chr_lengths)


## ----plot_map------------------------------------------------------------
# plot map
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])


## ----reorder_map---------------------------------------------------------
# reorder and flip some chromosomes
map2_chr_order <- c("Pet8", "Pet9", "Pet12_16", "Pet15", "Pet16_17", "Pet17")
flip_chrs <- c("Pet9")
map_list <- make_one_map(ann_pet_map, map1_max_chr_lengths = ann_chr_lengths, map2_chr_order = map2_chr_order, flip_chrs = flip_chrs)
plot_maps(map_df = map_list[[1]], map1_chrom_breaks = map_list[[2]], map2_chrom_breaks = map_list[[3]])


## ----find_params---------------------------------------------------------
# find best parameter combination
# run find_synteny_blocks with each parameter combination and collect summary statistics

parameter_data <- test_parameters(map_list, max_cluster_range_list = seq(1, 5, by = 0.5), max_nn_dist_list = seq(10, 50, by = 10))
plot_summary_stats(parameter_data[[1]], "composite")



## ----find_blocks---------------------------------------------------------
# find synteny blocks
synt_blocks <- find_synteny_blocks(map_list, max_cluster_range = 2, max_nn_dist = 10, plots = TRUE)


## ----find_blocks2--------------------------------------------------------

# print the contents of synt_blocks

lapply(synt_blocks, as_tibble)


## ----find_blocks3--------------------------------------------------------

# plot synteny block orientations
synt_blocks[[1]] %>%
    plot_maps(map_list[[2]], map_list[[3]], col = c("blue", "grey", "red", "black")[as.numeric(as.factor(.$orientation))], 
              main = "Synteny block orientation",  cex_val = 0.75)


