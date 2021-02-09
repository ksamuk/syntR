.onAttach <- function(...) {
  ver <- utils::packageVersion("syntR")
  cat(crayon::bold(crayon::cyan(paste0("Loaded syntR version ", ver, "!\n"))))
  cat(crayon::green(paste0("* Online documentation and vignettes at https://ksamuk.github.io/syntR/index.html\n")))
  cat(crayon::green(paste0("* Please cite the following paper if you use syntR in your work: \n")))
  cat(crayon::italic(crayon::green(paste0("* Kate L Ostevik, Kieran Samuk, Loren H Rieseberg, Ancestral Reconstruction of Karyotypes Reveals an Exceptional Rate of Nonrandom Chromosomal Evolution in Sunflower, Genetics, Volume 214, Issue 4, 1 April 2020, Pages 1031–1045, https://doi.org/10.1534/genetics.120.303026 \n"))))
}
