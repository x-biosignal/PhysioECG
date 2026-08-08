# Assess PPG signal quality

Computes population skewness and raw kurtosis, either per channel or in
non-overlapping windows.

## Usage

``` r
ppgQuality(x, window_sec = NULL, skew_threshold = -0.3, assay_name = NULL)
```

## Arguments

- x:

  A PPG `PhysioExperiment`.

- window_sec:

  Optional non-overlapping window duration in seconds.

- skew_threshold:

  Minimum skewness SQI accepted as usable.

- assay_name:

  Optional raw PPG assay name.

## Value

A data frame with channel, segment, start time, skewness SQI, kurtosis
SQI, and acceptance flag.

## References

Elgendi M (2016). Optimal signal quality index for photoplethysmogram
signals. *Bioengineering*, 3:21.
[doi:10.3390/bioengineering3040021](https://doi.org/10.3390/bioengineering3040021)

## Examples

``` r
simulation <- make_ppg(n_time = 2000)
ppgQuality(simulation$pe)
#>   channel segment start_sec    s_sqi    k_sqi accept
#> 1       1       1         0 1.998507 6.479899   TRUE
```
