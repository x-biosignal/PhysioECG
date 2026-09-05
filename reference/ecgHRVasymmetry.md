# Heart-Rate Asymmetry (Poincare Plot Asymmetry) Descriptors

Computes heart-rate asymmetry (HRA) descriptors from RR interval data.
The Poincare plot is asymmetric about the line of identity because
heart-rate *decelerations* (prolongations of RR, points above the line,
\\RR\_{n+1} \> RR_n\\) and *accelerations* (shortenings, points below
the line) contribute unequally to variability – a signature of vagal
activity and of the temporal irreversibility of heart-rate dynamics.
Decelerations typically dominate short-term variability while
accelerations dominate long-term variability (Piskorski & Guzik, 2011).

## Usage

``` r
ecgHRVasymmetry(rr)
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms`, and `time_sec`, as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

## Value

A data.frame with one row per channel and the columns:

- channel:

  Channel identifier.

- gi:

  Guzik's Index (\\ points from the line of identity divided by that of
  all points.

- si:

  Slope Index (\\ points divided by that of all points.

- ai:

  Area Index (\\ points divided by that of all points.

- pi:

  Porta's Index (\\ the number of points off the line of identity.

- c1d, c1a:

  Contributions of decelerations / accelerations to short-term HRV
  (`c1d + c1a = 1`).

- sd1d, sd1a:

  Short-term variance components of decelerations / accelerations (ms).

- c2d, c2a:

  Contributions of decelerations / accelerations to long-term HRV
  (`c2d + c2a = 1`).

- sd2d, sd2a:

  Long-term variance components (ms).

- cd, ca:

  Total contributions of decelerations / accelerations to HRV.

- sdnnd, sdnna:

  Total variance components (ms).

A value indicating asymmetry is `c1d > c1a` (decelerations dominate
short-term variability). Columns are `NA` for channels with fewer than
three RR intervals.

## References

Guzik, P., et al. (2006). "Heart rate asymmetry by Poincare plots of RR
intervals." *Biomedizinische Technik*, 51(4), 272–275.

Piskorski, J., & Guzik, P. (2011). "Asymmetric properties of long-term
and total heart rate variability." *Medical & Biological Engineering &
Computing*, 49(11), 1289–1297.
[doi:10.1007/s11517-011-0834-z](https://doi.org/10.1007/s11517-011-0834-z)

Porta, A., et al. (2008). "Temporal asymmetries of short-term heart
period variability are linked to autonomic regulation." *American
Journal of Physiology*, 295(2), R550–R557.

## See also

[`ecgHRVpoincare`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVpoincare.md)
for the symmetric Poincare descriptors (SD1/SD2),
[`ecgDFA`](https://x-biosignal.github.io/PhysioECG/reference/ecgDFA.md)
for detrended fluctuation analysis.
