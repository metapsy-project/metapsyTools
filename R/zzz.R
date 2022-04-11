.onAttach <- function(libname, pkgname) {
  packageStartupMessage(crayon::green("\u2713"),
                        crayon::cyan(" Loading "),
                        crayon::cyan$bold("{metapsyTools}"),
                        crayon::cyan(" 0.3.2 [BETA]. \n \u2192 For help, run "),
                        crayon::green("vignette('metapsyTools')"),
                        crayon::cyan(" in the R console."),
                        crayon::cyan("\n \u2192 Package documentation: "),
                        crayon::green("metapsytools.protectlab.org"),
                        crayon::cyan(".\n \u2192 Password: metapsyTools123!")) 
}



