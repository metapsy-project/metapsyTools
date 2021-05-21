# Simulate the number of treatment cycles and "excess treatments"

This function allows you to simulate the total number of treatment
cycles needed to achieve response in 99% of all patients in a
population, and the number of "excess" treatments due to non-response.

## Usage

``` r
simulateTreatmentCycles(response.rate, decay, max.cycles = 1000)
```

## Arguments

- response.rate:

  The response rate of the treatment(s). Can be either a single numeric
  value between 0.01 and 0.99, or a vector of single numeric values,
  indicating the response rate after the first, second, third, etc.
  treatment cycle.

- decay:

  The proportion (percentage) by which the response rate decreases after
  each unsuccessful treatment cycle. Must be a numeric value between
  -0.99 and 0.99. Negative values indicate an *increase* in response
  rates after each unsuccessful treatment cycle. A value of 0 indicates
  stable response rates in each cycle. If `response.rate` includes
  multiple values, the decay is applied after the last user-provided
  response rate.

- max.cycles:

  By default, the simulation is terminated if 99% cumulative response in
  the population has not been reached after 1000 cycles. This can be set
  to a lower number to model more realistic scenarios (e.g., no more
  treatment attempts after 10 cycles, so that `max.cycles=10`).

## Value

Returns an object of class `"simulateTreatmentCycles"`. This object
includes, among other things, a `data.frame` with the name `data`, in
which all simulation results are stored.

Other objects are the total number of treatments provided
(`total.treatments`), the total number of "excess" treatments
(`excess.treatments`), and the average number of treatment cycles per
person (`avg.no.treatments`), assuming a population of 100 patients. A
`plot` S3 method is also defined for the `simulateTreatmentCycles`
object.

## Details

In the “simple” scenario, the repeated treatment cycles can be modeled
as a Markov chain with absorbing state \\S_0\\ (patient responds), as
well as transition probabilities \\p_k=P(S_0│S_k)\\ and
\\1-p_k=P(S\_{k+1}\|S_k)\\. Here, \\k\\ is the current treatment cycle,
and \\p_k\\ is the probability of responding to the \\k\\th treatment,
which is held constant across all treatment cycles. In this scenario,
the formula to calculate the cumulative response \\C_n\\ after \\n\\
treatment cycles reduces to:

\\C_n=1-(1-p)^n\\

The number of “excess” treatments \\N^E\\ can be obtained using this
formula (assuming that 100 patients are treated initially):

\\N^E\approx\frac{100×(1-(1-p)^n)}{p}-100\\

If we additionally consider “decay” of the treatment response (e.g.,
response rates decrease by 10% with each additional treatment attempt),
we obtain the following formula for the cumulative response:

\\C_n = \sum\_{i=1}^{n} \left( \prod\_{j=1}^{i-1} \left( 1 -
p(1-d)^{j-1} \right) \right) p(1-d)^{i-1}\\

where \\C_n\\ is the cumulative treatment response, \\p\\ is the
response rate of treatments, \\d\\ encodes the proportional “decay” with
each treatment cycle, and \\n\\ is the total number of treatment cycles.

For more details on the metapsyTools package, see the [Get
Started](https://tools.metapsy.org/articles/metapsytools) vignette.

## Author

Mathias Harrer <mathias.h.harrer@gmail.com>, Pim Cuijpers
<p.cuijpers@vu.nl>, Toshi Furukawa

## Examples

``` r
if (FALSE) { # \dontrun{
# "Simple" scenario: 50% response rate at each cycle
res <- simulateTreatmentCycles(response.rate = 0.5, decay = 0)
res; plot(res)

# Define manual response rates. Since decay=0, rates stay at 0.4 afterwards
res <- simulateTreatmentCycles(response.rate = c(0.57, 0.4), decay = 0)
res; plot(res)

# Define manual response rates up to the 4th cycles. Afterwards, response decays by 10%
res <- simulateTreatmentCycles(response.rate = c(0.56, 0.3, 0.28, 0.4), decay = .1)
res; plot(res)

# Run the simulation for a fixed number of cycles. Here: 10
res <- simulateTreatmentCycles(response.rate = c(0.56, 0.3, 0.28, 0.1), decay = .1, max.cycles = 10)
res; plot(res)

# Go crazy
res <- simulateTreatmentCycles(response.rate = c(.2, .8, .2, .8, .2, .8, .2), decay = -.2)
res; plot(res)
} # }
```
