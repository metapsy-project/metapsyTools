# Calculate the log ratio of means (ROM) using within-group change data

Computes the natural log-transformed ratio of means (log ROM) and its
standard error using the delta method, based on pre-post change scores
rather than post-treatment means. This is an internal helper called
row-wise by [`calculateEffectSizes`](calculateEffectSizes.md) and is not
intended to be used directly.

## Usage

``` r
rom.change.m.sd(x, ...)
```

## Arguments

- x:

  A `data.frame` in which each row represents one trial arm comparison.
  Must contain the columns listed under `...`.

- ...:

  The following columns are required and consumed from `x`:

  `mean_change_arm1`

  :   Mean change from baseline in arm 1 (treatment).

  `mean_change_arm2`

  :   Mean change from baseline in arm 2 (control/comparator).

  `sd_change_arm1`

  :   Standard deviation of the change score in arm 1.

  `sd_change_arm2`

  :   Standard deviation of the change score in arm 2.

  `n_arm1`

  :   Sample size of arm 1. Must be \> 0.

  `n_arm2`

  :   Sample size of arm 2. Must be \> 0.

  Additional columns in `x` are silently ignored.

## Value

A `data.frame` with the same number of rows as `x` and two numeric
columns:

- `es`:

  Log ratio of change-score means, \\\ln(\bar{x}\_{\Delta,1} /
  \bar{x}\_{\Delta,2})\\.

- `se`:

  Standard error of `es`, derived via the delta method.

Rows are set to `NA` for both columns when any required input is `NA`,
when either change-score mean equals zero (log ROM is undefined), when
either sample size is non-positive, or when the computed variance is
non-finite or not strictly positive.

## Details

The estimator is identical in form to that used in
[`rom.m.sd`](rom.m.sd.md), but the inputs are mean change scores and
their standard deviations: \$\$\ln(\text{ROM}) =
\ln\\\left(\frac{\bar{x}\_{\Delta,1}}{\bar{x}\_{\Delta,2}}\right)\$\$

Its sampling variance is approximated via the delta method (Hedges et
al., 1999): \$\$v = \frac{s\_{\Delta,1}^2}{n_1
\\\bar{x}\_{\Delta,1}^2} + \frac{s\_{\Delta,2}^2}{n_2\\
\bar{x}\_{\Delta,2}^2}\$\$

The standard error returned in `se` is \\\sqrt{v}\\.

A positive log ROM indicates a larger mean change in arm 1 than in arm
2. The measure is undefined when either change-score mean equals zero,
so such rows are returned as `NA`. Note that a zero mean change score is
substantively meaningful (no average change from baseline) and
relatively common, so this limitation should be considered when choosing
between effect size metrics.

## References

Hedges, L. V., Gurevitch, J., & Curtis, P. S. (1999). The meta-analysis
of response ratios in experimental ecology. *Ecology, 80*(4), 1150–1156.
[doi:10.1890/0012-9658(1999)080\[1150:TMAORR\]2.0.CO;2](https://doi.org/10.1890/0012-9658%281999%29080%5B1150%3ATMAORR%5D2.0.CO%3B2)

Lajeunesse, M. J. (2011). On the meta-analysis of response ratios for
studies with correlated and multi-group designs. *Ecology, 92*(11),
2049–2055. [doi:10.1890/11-0423.1](https://doi.org/10.1890/11-0423.1)

## See also

[`rom.m.sd`](rom.m.sd.md) for the post-treatment mean variant;
[`calculateEffectSizes`](calculateEffectSizes.md) for the calling
function.
