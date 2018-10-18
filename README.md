
<!-- README.md is generated from README.Rmd. Please edit that file -->
syntR <img src="inst/figures/logo.png" align="right" width="120" height="135" />
================================================================================

syntR is an R package for the reproducible identification of synteny blocks/chromosomal rearrangments via comparison of two genetic maps. syntR implements an error-aware clustering algorithm specifically designed for the highly linear structure of comparative genetic map data. syntR can be used to identify synteny blocks using any type of ordered and/or aligned genetic markers. 

***Note: syntR is currently in active development.***

Installation
------------

You can install syntR from github with:

``` r
# install.packages("devtools")
devtools::install_github("ksamuk/syntR")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```

Authors
------------
[Katherine Ostevik](http://www.kateostevik.com/) and [Kieran Samuk](https://ksamuk.github.io/).

See Also
-------
[GRIMM](http://grimm.ucsd.edu/GRIMM/) - A tool for analyzing rearrangements in pairs of genomes, including unichromosomal and multichromosomal genomes, and signed and unsigned data. 

[Flagel et al. 2018](https://www.biorxiv.org/content/early/2018/05/26/330159) - An example of a more formal model-based approach to a similar problem. 
