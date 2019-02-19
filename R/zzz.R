.onAttach <- function(...) {
  ver <- utils::packageVersion("syntR")
  packageStartupMessage("This is syntR version ", ver)
  packageStartupMessage("* Online documentation and vignettes at https://ksamuk.github.io/syntR/index.html")
  packageStartupMessage("* Please cite the following paper if you use syntR in your work: ")
  packageStartupMessage("* [citation forthcoming] ")
}
