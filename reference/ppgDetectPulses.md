# Detect systolic peaks in PPG

Implements the Elgendi two event-related moving-average (TERMA)
detector. The shipped windowed-sinc FIR is used in place of the original
Butterworth filter. Its linear-phase group delay is explicitly undone by
snapping each event-block candidate to the corresponding raw-signal
maximum.

## Usage

``` r
ppgDetectPulses(
  x,
  w1_ms = 111,
  w2_ms = 667,
  beta = 0.02,
  low = 0.5,
  high = 8,
  order = 65,
  refractory_ms = 300,
  assay_name = NULL
)
```

## Arguments

- x:

  A PPG `PhysioExperiment`.

- w1_ms:

  Systolic event window in milliseconds.

- w2_ms:

  Beat window in milliseconds.

- beta:

  Moving-average threshold offset.

- low, high:

  PPG band-pass cutoffs in Hz.

- order:

  FIR order; even values are promoted to the next odd value.

- refractory_ms:

  Minimum systolic-peak separation in milliseconds.

- assay_name:

  Optional raw PPG assay name.

## Value

A data frame with `channel`, `sample`, `time_sec`, and raw `amplitude`.

## Details

This implementation uses the package's pure-R windowed-sinc FIR rather
than the Butterworth filter in the original method. A strong late
dicrotic wave may require increasing `refractory_ms` so that it is not
counted as a second systolic pulse.

## References

Elgendi M, Norton I, Brearley M, Abbott D, Schuurmans D (2013). Systolic
peak detection in acceleration photoplethysmograms measured from
emergency responders in tropical conditions. *PLoS ONE*, 8:e76585.
[doi:10.1371/journal.pone.0076585](https://doi.org/10.1371/journal.pone.0076585)

## Examples

``` r
set.seed(1)
simulation <- make_ppg(n_time = 3000)
head(ppgDetectPulses(simulation$pe))
#>   channel sample time_sec amplitude
#> 1       1    122    0.968  2.996094
#> 2       1    226    1.800  3.005296
#> 3       1    331    2.640  2.993458
#> 4       1    435    3.472  3.011860
#> 5       1    539    4.304  3.012044
#> 6       1    643    5.136  3.014261
```
