# Plot HRV Power Spectral Density with Frequency Bands

Plots the RR-interval power spectrum with the VLF, LF and HF bands
shaded and their integrated powers annotated. The PSD and band
integration reuse the same internals as
[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md),
so the shaded band areas equal the `ecgHRVfreq` band powers. The band
powers are attached as the `"band_power"` attribute.

## Usage

``` r
plotHRVpsd(
  rr,
  method = c("welch", "ar", "lomb"),
  detrend = FALSE,
  channel = NULL,
  vlf_band = c(0.003, 0.04),
  lf_band = c(0.04, 0.15),
  hf_band = c(0.15, 0.4)
)
```

## Arguments

- rr:

  An rr data.frame or an object from
  [`ecgProcess`](https://x-biosignal.github.io/PhysioECG/reference/ecgProcess.md).

- method:

  Spectral method ("welch" default, "ar" or "lomb"), matching
  [`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md).

- detrend:

  Smoothness-priors detrending, matching
  [`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md).

- channel:

  Channel to plot (default: first).

- vlf_band, lf_band, hf_band:

  Band edges (Hz).

## Value

A ggplot object with a `"band_power"` attribute.

## References

Task Force (1996). Heart rate variability. *Circulation*, 93(5),
1043-1065.

## See also

[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md)

## Examples

``` r
pe <- make_ecg(40000, sr = 500)
rr <- ecgRRintervals(pe, ecgDetectRpeaks(pe))
plotHRVpsd(rr)
```
