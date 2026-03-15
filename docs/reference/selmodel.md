# Get or set selection model steps

Generic getter and setter for the selection model steps
(`selmodelSteps`) component of a `runMetaAnalysis` object (after
`correctPublicationBias`). See
[`Replacement functions`](Replacement-functions.md).

## Usage

``` r
selmodel(x, ...)

selmodel(x) <- value
```

## Arguments

- x:

  Object of class `runMetaAnalysis`.

- ...:

  Optional arguments (ignored).

- value:

  New value for `selmodelSteps`.

## Value

For `selmodel(x)`, the current `selmodelSteps` value. For
`selmodel(x) <- value`, the updated object (invisibly).
