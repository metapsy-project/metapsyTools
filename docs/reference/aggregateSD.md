# Aggregate standard deviations of multiple trials

This function allows you to aggregate standard deviations of outcome
scores across multiple trials.

## Usage

``` r
aggregateSD(sd, n, ...)
```

## Arguments

- sd:

  A numeric vector of standard deviations.

- n:

  A numeric vector of sample sizes corresponding to each standard
  deviation.

- ...:

  Additional arguments (not used).

## Value

The pooled standard deviation across trials.

## Details

Values are aggregated assuming homogeneous variances across studies
(i.e., all trials share the same population standard deviation). The
following formula is used:

\\s\_{\text{pooled}} = \sqrt{ \frac{ \sum (n_i - 1) \cdot s_i^2 }{ \sum
(n_i - 1) } }\\

where \\s_i\\ and \\n_i\\ are the standard deviation and sample size in
some study \\i\\. Please note that aggregation is only possible when all
SD values come from the same scale or instrument.

## Examples

``` r
if (FALSE) { # \dontrun{
sd <- c(8.77, 9.01, 10.11)
n <- c(120, 15, 80)
aggregateSD(sd, n)
} # }
```
