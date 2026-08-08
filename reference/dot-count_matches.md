# Count template matches for sample entropy

Counts the number of template vector pairs of length `dim` that are
within tolerance `r` (Chebyshev distance). Uses a vectorised approach
with an embedding matrix to avoid nested R loops for better performance
on long RR series.

## Usage

``` r
.count_matches(x, dim, r)
```

## Arguments

- x:

  Numeric vector (time series).

- dim:

  Embedding dimension.

- r:

  Tolerance.

## Value

Integer count of matching template pairs.
