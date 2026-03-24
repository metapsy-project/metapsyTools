# Calculate NNTs

Calculate NNTs (extracted from `dmetar`)

## Usage

``` r
metapsyNNT(d, CER, event.e, n.e, event.c, n.c, RR, names, method)
```

## Arguments

- d:

  A single numeric or concatenated vector of numerics representing the
  effect size expressed as Cohen's \\d\\ or Hedges' \\g\\. If this is
  the only parameter specified in the function, the method by Kraemer
  and Kupfer is used automatically to calculate \\NNT\\s.

- CER:

  The control group event ratio. Furukawa's method (Furukawa &
  Leucht, 2011) to calculate `NNT`s from `d` requires that the assumed
  response ("event") ratio in the control group
  (\\\frac{n\_{responders}}{N\_{total}}\\) is specified. The CER can
  assume values from 0 to 1. If a value is specified for `CER`,
  Furukawa's method is used automatically. Argument `method` has to be
  set to `"KraemerKupfer"` to override this.

- event.e:

  Single number or numeric vector. The number of (favourable) events in
  the experimental group.

- n.e:

  Single number or numeric vector. The number participants in the
  experimental group.

- event.c:

  Single number or numeric vector. The number of (favourable) events in
  the control group.

- n.c:

  Single number or numeric vector. The number of participants in the
  control group.

- RR:

  Optional. A single numeric or vector of numeric risk ratios. When
  supplied, NNT is calculated from `RR` and `CER` via \\NNT = 1 / \|CER
  \* (RR - 1)\|\\. This requires `CER` to be provided as a numeric value
  between 0 and 1.

- names:

  Optional. Character vector of equal length as the vector supplied to
  `d` or `event.e` containing study/effect size labels.

- method:

  The method to be used to calculate the NNT from `d`. Either
  `"KraemerKupfer"` for the method proposed by Kraemer and Kupfer (2006)
  or `"Furukawa"` for the Furukawa method (Furukawa & Leucht, 2011).
  Please note that the Furukawa's method can only be used when `CER` is
  specified.
