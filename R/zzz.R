.onAttach <- function(libname, pkgname) {
  packageStartupMessage(crayon::green("\u2713"),
                        crayon::cyan(" Loading "),
                        crayon::cyan$bold("{metapsyTools}"),
                        crayon::cyan(" 0.2.1 [BETA]. \n \u2192 For help, run "),
                        crayon::green("vignette('metapsyTools')"),
                        crayon::cyan(" in the R console."))
}



