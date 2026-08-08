# Plot Annotated ECG Waveform

Plots the ECG waveform with P, QRS and T fiducials from
[`ecgDelineate`](https://x-biosignal.github.io/PhysioECG/reference/ecgDelineate.md).
The detected R-peak sample indices are attached as the `"peaks"`
attribute (a beat-editing hook); the delineation is attached as
`"delineation"`.

## Usage

``` r
plotEcgWave(x, peaks, channel = 1, window_sec = 10, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- peaks:

  R-peak table from
  [`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md).

- channel:

  Channel to plot (default: 1).

- window_sec:

  Seconds of signal to display from the start (default 10; NULL for the
  whole recording).

- assay_name:

  Input assay name (default: first assay).

## Value

A ggplot object with `"peaks"` and `"delineation"` attributes.

## References

Pan, J. & Tompkins, W.J. (1985). A real-time QRS detection algorithm.
*IEEE TBME*, 32(3), 230-236.

## See also

[`ecgDelineate`](https://x-biosignal.github.io/PhysioECG/reference/ecgDelineate.md),
[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)

## Examples

``` r
pe <- make_ecg(5000, sr = 500)
pk <- ecgDetectRpeaks(pe)
plotEcgWave(pe, pk)
```
