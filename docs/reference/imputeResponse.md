# Impute response rates based on continuous outcome data

This function allows to impute response rates based on the post-test
means, standard deviations and sample size using the
normal-approximation method by Furukawa et al. (2005). Response rates
can be imputed using the 50% symptom decrease threshold, user-defined
cut-off values, or the reliable change index (Jacobson & Truax, 1991).

## Usage

``` r
imputeResponse(m.trt.pre, m.trt.post, sd.trt.post, n.trt, 
               m.ctr.pre, m.ctr.post, sd.ctr.post, n.ctr,
               sd.trt.pre, sd.ctr.pre, 
               rho = 0.8, cutoff = NULL, lower.is.better = TRUE,
               method = c("cutpoint", "rci"))
```

## Arguments

- m.trt.pre:

  Pre-test mean in the (treatment) arm.

- m.trt.post:

  Post-test mean in the (treatment) arm.

- sd.trt.post:

  Post-test standard deviation in the (treatment) arm.

- n.trt:

  Sample size in the (treatment) arm.

- m.ctr.pre:

  *Optional*. Pre-test mean in the control arm.

- m.ctr.post:

  *Optional*. Post-test mean in the control arm.

- sd.ctr.post:

  *Optional*. Post-test standard deviation in the control arm.

- n.ctr:

  *Optional*. Sample size in the control arm.

- sd.trt.pre:

  Pre-test standard deviation in the (treatment) group. Required if
  `method = "rci"`.

- sd.ctr.pre:

  *Optional*. Pre-test standard deviation in the control group. Required
  if `method = "rci"`.

- rho:

  Reliability coefficient, used to calculate the reliable change when
  `method = "rci"`. Set to 0.8 by default.

- cutoff:

  User defined cut-off score for response. `NULL` by default, which
  means that 50% symptom decrease is used as the response criterion.

- lower.is.better:

  Do lower values indicate better outcomes (e.g. less depression)?
  `TRUE` by default.

- method:

  Method to define response. Either `"cutpoint"` (default; uses 50%
  symptom decrease as criterion, or a user-specified score in `cutoff`)
  or `"rci"`, which means that the reliable change index (Jacobson &
  Truax, 1991) is used.

## Value

If only values for the treatment arm are specified, `imputeResponse`
returns a `data.frame` with two columns:

- `trtResponder`: Number of responders

- `nTrt`: Total sample size of the treatment arm

If values are also specified for the control group, six additional
columns are returned:

- `ctrResponder`: Number of responders in the control group

- `nCtr`: Total sample size of the treatment arm

- `logRR`: Log-risk ratio of the difference in response rates in both
  arms

- `seLogRR`: Standard error of the log-risk ratio

- `logOR`: Log-odds ratio of the difference in response rates in both
  arms

- `seLogOR`: Standard error of the log-odds ratio

## Details

For more details see the [Get
Started](https://tools.metapsy.org/articles/metapsytools) vignette.

## References

Furukawa, T. A., Cipriani, A., Barbui, C., Brambilla, P., & Watanabe, N.
(2005). Imputing response rates from means and standard deviations in
meta-analyses. *International Clinical Psychopharmacology, 20*(1),
49-52.

Jacobson, N. S., & Truax, P. (1991). Clinical significance: a
statistical approach to defining meaningful change in psychotherapy
research. *Journal of Consulting and Clinical Psychology, 59*(1), 12–19.
https://doi.org/10.1037//0022-006x.59.1.12

## See also

[`calculateEffectSizes`](calculateEffectSizes.md)

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{

# Calculate response using 50% method for one group.
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100)

# Calculate response for two groups and calculate effect sizes
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120)

# Calculate using user-defined response threshold
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120,
               cutoff = 15)

# Assuming higher outcomes are better...
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120,
               cutoff = 15, lower.is.better = FALSE)

# Using RCI, with an assumed reliability of 0.88
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120,
               method = "rci", rho = 0.88,
               sd.trt.pre = 8)

# Using a different pre-test SD in the control
imputeResponse(m.trt.pre = 20, m.trt.post = 11,
               sd.trt.post = 7, n.trt = 100,
               m.ctr.pre = 20, m.ctr.post = 18,
               sd.ctr.post = 6, n.ctr = 120,
               method = "rci", rho = 0.88,
               sd.trt.pre = 8, sd.ctr.pre = 9)

# Calculating many results at the same time
set.seed(123)
imputeResponse(m.trt.pre = runif(10, 10, 30), m.trt.post = runif(10, 8, 14),
               sd.trt.post = runif(10, 5, 8), n.trt = rpois(10, 150),
               m.ctr.pre = runif(10, 10, 30), m.ctr.post = runif(10, 13, 20),
               sd.ctr.post = runif(10, 5, 8), n.ctr = rpois(10, 150),
               method = "rci", rho = round(runif(10, 70, 90))/100,
               sd.trt.pre = runif(10, 7, 9))

set.seed(123)
imputeResponse(m.trt.pre = runif(10, 10, 30), m.trt.post = runif(10, 8, 14),
               sd.trt.post = runif(10, 5, 8), n.trt = rpois(10, 150),
               m.ctr.pre = runif(10, 10, 30), m.ctr.post = runif(10, 13, 20),
               sd.ctr.post = runif(10, 5, 8), n.ctr = rpois(10, 150))

set.seed(123)
imputeResponse(m.trt.pre = runif(10, 10, 30), m.trt.post = runif(10, 8, 14),
               sd.trt.post = runif(10, 5, 8), n.trt = rpois(10, 150),
               cutoff = 11:20)
} # }
```
