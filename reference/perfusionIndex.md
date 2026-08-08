# Calculate the PPG perfusion index

Perfusion index is 100 times the pulsatile systolic-to-foot amplitude
divided by the raw-signal DC level.

## Usage

``` r
perfusionIndex(x, peaks = NULL, assay_name = NULL, agg = c("median", "mean"))
```

## Arguments

- x:

  A raw, DC-preserving PPG `PhysioExperiment`.

- peaks:

  Optional systolic peak table from
  [`ppgDetectPulses()`](https://x-biosignal.github.io/PhysioECG/reference/ppgDetectPulses.md).

- assay_name:

  Optional raw PPG assay name.

- agg:

  Aggregation of per-pulse AC amplitudes.

## Value

A per-channel data frame containing AC, DC, perfusion index, and pulse
count.

## References

Charlton PH, Kyriacou PA, Mant J, et al. (2022). Wearable
photoplethysmography for cardiovascular monitoring. *Proceedings of the
IEEE*, 110:355-381.
[doi:10.1109/JPROC.2022.3149785](https://doi.org/10.1109/JPROC.2022.3149785)

## Examples

``` r
simulation <- make_ppg(n_time = 5000)
perfusionIndex(simulation$pe)
#>   channel       ac       dc   pi_pct n_pulses
#> 1       1 1.097512 2.148277 51.08803       47
```
