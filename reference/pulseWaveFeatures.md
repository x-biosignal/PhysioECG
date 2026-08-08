# Extract PPG pulse-wave features

Measures foot-to-systolic rise time, pulse amplitude, dicrotic and
diastolic fiducials, reflection and augmentation indices, pulse width,
and optional ECG-to-PPG pulse transit time and pulse-wave velocity.

## Usage

``` r
pulseWaveFeatures(
  x,
  peaks = NULL,
  ecg_peaks = NULL,
  distance_m = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A raw PPG `PhysioExperiment`.

- peaks:

  Optional systolic peak table from
  [`ppgDetectPulses()`](https://x-biosignal.github.io/PhysioECG/reference/ppgDetectPulses.md).

- ecg_peaks:

  Optional ECG R-peak table containing `sample` and optionally
  `channel`.

- distance_m:

  Optional ECG-to-PPG path length in metres for PWV.

- assay_name:

  Optional raw PPG assay name.

## Value

A data frame with one row per pulse and waveform/PTT/PWV features.

## References

Elgendi M (2012). On the analysis of fingertip photoplethysmogram
signals. *Current Cardiology Reviews*, 8:14-25.
[doi:10.2174/157340312801215782](https://doi.org/10.2174/157340312801215782)

## Examples

``` r
simulation <- make_ppg(n_time = 5000)
head(pulseWaveFeatures(simulation$pe))
#>   channel beat foot_sample systolic_sample systolic_time_sec rise_time_ms
#> 1       1    1         103             123             0.976          160
#> 2       1    2         208             227             1.808          152
#> 3       1    3         312             330             2.632          144
#> 4       1    4         416             436             3.480          160
#> 5       1    5         520             540             4.312          160
#> 6       1    6         624             644             5.144          160
#>   pulse_amplitude dicrotic_notch_sample diastolic_sample reflection_index
#> 1        1.093694                   133              147         47.08573
#> 2        1.115187                   238              254         46.19351
#> 3        1.077151                   331              332         97.32652
#> 4        1.113700                   447              459         46.70862
#> 5        1.074340                   550              563         46.67595
#> 6        1.105202                   655              667         45.27318
#>   augmentation_index pulse_width_ms ptt_ms pwv_m_s
#> 1          52.914269            840     NA      NA
#> 2          53.806492            832     NA      NA
#> 3           2.673476            832     NA      NA
#> 4          53.291382            832     NA      NA
#> 5          53.324054            832     NA      NA
#> 6          54.726817            832     NA      NA
```
