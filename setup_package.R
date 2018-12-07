# setup pacakge

library("devtools")
library("usethis")
library("roxygen2")

# git/github
# set up a new github personal access token
use_git()
browse_github_pat() # interactive
edit_r_environ() # interactive
use_github()

# dependencies
use_package("dplyr")
use_package("tidyr")
use_package("reshape2")
use_package("intervals")
use_package("viridisLite")

# data
ann_pet_map <- read.table("data/tutorial_data_ann_pet_comparison.txt", h = T)
ann_chr_lengths <- read.table("data/tutorial_data_ann_chr_max.txt", h = T)
use_data(ann_pet_map)
use_data(ann_chr_length)
#use_data(ann_chr_length)

# documentation
use_readme_rmd()
use_package_doc()
document()

# vignenette
use_vignette("syntR_tutorial")
build_vignettes()
browseVignettes("syntR")

# load package
load_all()
