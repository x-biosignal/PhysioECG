# Accept/reject ECG channels by bSQI

Thin wrapper over
[`ecgBeatSQI`](https://x-biosignal.github.io/PhysioECG/reference/ecgBeatSQI.md)
adding a logical `accept` column for channels whose bSQI meets a
threshold.

## Usage

``` r
ecgBeatSQIgate(x, threshold = 0.8, ...)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- threshold:

  Minimum bSQI to accept a channel (default 0.8).

- ...:

  Passed to
  [`ecgBeatSQI`](https://x-biosignal.github.io/PhysioECG/reference/ecgBeatSQI.md).

## Value

The
[`ecgBeatSQI`](https://x-biosignal.github.io/PhysioECG/reference/ecgBeatSQI.md)
data.frame with an added logical `accept` column.

## See also

[`ecgBeatSQI`](https://x-biosignal.github.io/PhysioECG/reference/ecgBeatSQI.md)

## Examples

``` r
pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
ecgBeatSQIgate(pe, threshold = 0.8)
#>   channel n_a n_b n_matched bSQI accept
#> 1       1  11  11        11    1   TRUE
```
