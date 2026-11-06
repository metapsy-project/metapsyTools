# Add effect size columns needed for BLUP calculations to the model data

`metabin` models do not store `.TE` and `.seTE` in `$data`, so these are
always taken from the model itself.

## Usage

``` r
blupData(M)
```
