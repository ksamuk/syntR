# setup pacakge

library("devtools")
library("usethis")
library("roxygen2")

# git/github
use_git()
# set up a new github personal access token
# after the 'cat', restart Rstudio
browse_github_pat() # interactive
edit_r_environ() # interactive
use_github()

# dependencies
use_package("dplyr")
use_package("tidyr")
use_package("reshape2")
use_package("MASS")
use_package("intervals")
use_package("viridisLite")

# data
ann_pet_map <- read.table("data/ann_pet_sample_comparison.txt", h = T)
ann_chr_length <- read.table("data/Ann_chr_max.txt", h = F) # needs chrosome names?
use_data(ann_pet_map)
use_data(ann_chr_length)
#use_data(ann_chr_length)

# documentation
document()

# load package
load_all()
