.onAttach <- function(libname, pkgname) {
  packageStartupMessage(crayon::green("\u2713"),
                        crayon::cyan(" Loading "),
                        crayon::cyan$bold("{metapsyTools}"),
                        crayon::cyan(" 1.0.4 [BETA]."),
                        crayon::cyan("\n \u2192 Package documentation: "),
                        crayon::green("tools.metapsy.org")) 
}



