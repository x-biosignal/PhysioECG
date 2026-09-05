# Cardiac Autonomic Indices (Toichi CSI / CVI) from the Poincare Plot

Computes the cardiac sympathetic and vagal indices (Toichi et al., 1997)
from RR interval data. The indices are read off the Poincare plot's
\\4\times\\SD bounding box: the long axis \\L = 4\\\mathrm{SD2}\\ and
the short axis \\T = 4\\\mathrm{SD1}\\.

- `csi`:

  Cardiac Sympathetic Index, \\L/T\\ (= SD2/SD1); rises with sympathetic
  predominance.

- `cvi`:

  Cardiac Vagal Index, \\\log\_{10}(L\times T)\\; rises with
  parasympathetic (vagal) activity.

- `csi_modified`:

  Modified CSI, \\L^2/T\\, a more sensitive sympathetic marker (Jeppesen
  et al., 2014) used in seizure detection.

Here SD1 is the universal beat-to-beat term (\\\mathrm{SDSD}/\sqrt2\\)
and SD2 is the geometric (paired-projection) long-term term, matching
the convention in which the Toichi indices are defined (as in
NeuroKit2's `hrv_nonlinear`); this SD2 differs by a small \\O(1/n)\\
amount from the analytical closed-form SD2 reported by
[`ecgHRVpoincare`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVpoincare.md).

## Usage

``` r
ecgHRVautonomic(rr)
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms`, and `time_sec`, as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

## Value

A data.frame with one row per channel and the columns `channel`, `csi`,
`cvi`, and `csi_modified`. Columns are `NA` for channels with fewer than
three RR intervals.

## References

Toichi, M., Sugiura, T., Murai, T., & Sengoku, A. (1997). "A new method
of assessing cardiac autonomic function and its comparison with spectral
analysis and coefficient of variation of R-R interval." *Journal of the
Autonomic Nervous System*, 62(1-2), 79–84.

Jeppesen, J., et al. (2014). "Using Lorenz plot and Cardiac Sympathetic
Index of heart rate variability for detecting seizures for patients with
epilepsy." *IEEE EMBC*, 4563–4566.

## See also

[`ecgHRVpoincare`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVpoincare.md)
for SD1/SD2,
[`ecgHRVasymmetry`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVasymmetry.md)
for heart-rate asymmetry.
