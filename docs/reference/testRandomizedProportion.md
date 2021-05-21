# Test if the randomized proportion differs from the original allocation ratio

This function returns a one-sample proportion test that allows to
examine if the included sample in a study deviated significantly form
the intended allocation ratio.

## Usage

``` r
testRandomizedProportion(
  n_arm1,
  n_arm2,
  study,
  data = NULL,
  ratio = "1:1",
  ...
)
```

## Arguments

- n_arm1:

  A vector in which the sample size of the first trial arm is stored.

- n_arm2:

  A vector in which the sample size of the second trial arm is stored.

- study:

  An optional vector with study labels.

- data:

  An optional `data.frame` in which the sample size data is included.

- ratio:

  A single string or character vector of the expected allocation ratio.
  This ratio should be separated by a colon, e.g. `"1:3`, where the left
  side indicates the first trial arm (`n_arm1`) and the right side
  indicates the second trial arm (`n_arm2`). Defaults to `"1:1`.

- ...:

  Additional arguments.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data that follows the Metapsy data standard
data("depressionPsyCtr")

# This is an unexported function, so we need the metapsyTools::: prefix
# Test if all studies follow the 1:1 allocation ratio
metapsyTools:::testRandomizedProportion(
  n_arm1, n_arm2, data = depressionPsyCtr, study = study) 

# Provide a comparison-specific allocation ratio
ratio <- c("1:1", "3:1", "4:5")
metapsyTools:::testRandomizedProportion(
  n_arm1, n_arm2, data = depressionPsyCtr[1:3,], 
  study = study, ratio = ratio) 
} # }
```
