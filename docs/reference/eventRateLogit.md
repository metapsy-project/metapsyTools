# Logit-transformed event rates and their standard errors

Event rates of 0 and 1 have an undefined logit. In these cases, a
continuity correction of 0.5 is added to the number of events and
non-events (as in
[`meta::metaprop`](https://rdrr.io/pkg/meta/man/metaprop.html) or
[`metafor::escalc`](https://wviechtb.github.io/metafor/reference/escalc.html)).

## Usage

``` r
eventRateLogit(event, n)
```
