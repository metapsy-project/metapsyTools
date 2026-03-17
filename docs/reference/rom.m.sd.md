# Calculate the log ratio of means (ROM) using means and standard deviations

Computes the natural log-transformed ratio of means (log ROM, also
called the log response ratio) and its standard error using the delta
method. This is an internal helper called row-wise by
[`calculateEffectSizes`](calculateEffectSizes.md) and is not intended to
be used directly.

## Usage

``` r
rom.m.sd(x, ...)
```

## Arguments

- x:

  A `data.frame` in which each row represents one trial arm comparison.
  Must contain the columns listed under `...`.

- ...:

  The following columns are required and consumed from `x`:

  `mean_arm1`

  :   Mean of the outcome in arm 1 (treatment).

  `mean_arm2`

  :   Mean of the outcome in arm 2 (control/comparator).

  `sd_arm1`

  :   Standard deviation of the outcome in arm 1.

  `sd_arm2`

  :   Standard deviation of the outcome in arm 2.

  `n_arm1`

  :   Sample size of arm 1. Must be \> 0.

  `n_arm2`

  :   Sample size of arm 2. Must be \> 0.

  Additional columns in `x` are silently ignored.

## Value

A `data.frame` with the same number of rows as `x` and two numeric
columns:

- `es`:

  Log ratio of means, \\\ln(\bar{x}\_1 / \bar{x}\_2)\\.

- `se`:

  Standard error of `es`, derived via the delta method.

Rows are set to `NA` for both columns when any required input is `NA`,
when either mean equals zero (log ROM is undefined), when either sample
size is non-positive, or when the computed variance is non-finite or not
strictly positive.

## Details

The log ROM effect size is defined as: \$\$\ln(\text{ROM}) =
\ln\\\left(\frac{\bar{x}\_1}{\bar{x}\_2}\right)\$\$

Its sampling variance is approximated via the delta method (Hedges et
al., 1999): \$\$v = \frac{s_1^2}{n_1 \bar{x}\_1^2} + \frac{s_2^2}{n_2
\bar{x}\_2^2}\$\$

The standard error returned in `se` is \\\sqrt{v}\\.

A positive log ROM indicates that arm 1 has a higher mean than arm 2.
The measure is undefined when either group mean equals zero, which is
why such rows are returned as `NA`.

## References

Hedges, L. V., Gurevitch, J., & Curtis, P. S. (1999). The meta-analysis
of response ratios in experimental ecology. *Ecology, 80*(4), 1150–1156.
[doi:10.1890/0012-9658(1999)080\[1150:TMAORR\]2.0.CO;2](https://doi.org/10.1890/0012-9658%281999%29080%5B1150%3ATMAORR%5D2.0.CO%3B2)

Lajeunesse, M. J. (2011). On the meta-analysis of response ratios for
studies with correlated and multi-group designs. *Ecology, 92*(11),
2049–2055. [doi:10.1890/11-0423.1](https://doi.org/10.1890/11-0423.1)

## See also

[`calculateEffectSizes`](calculateEffectSizes.md)
