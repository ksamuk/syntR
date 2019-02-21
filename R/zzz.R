.onAttach <- function(...) {
  ver <- utils::packageVersion("syntR")
  cat(crayon::bold(crayon::cyan(paste0("Loaded syntR version ", ver, "!\n"))))
  cat(crayon::green(paste0("* Online documentation and vignettes at https://ksamuk.github.io/syntR/index.html\n")))
  cat(crayon::green(paste0("* Please cite the following paper if you use syntR in your work: \n")))
  cat(crayon::italic(crayon::green(paste0("* [citation forthcoming] \n"))))
}
