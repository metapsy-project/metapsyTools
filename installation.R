# # # # # # # # #
# INSTALLATION  #
# # # # # # # # #

# Run this code:

if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("MathiasHarrer/metapsyTools",
                         build_vignettes = TRUE,
                         auth_token = "0acdd3b1fb4e36cdc7efa5f59da2918d121ea051")
