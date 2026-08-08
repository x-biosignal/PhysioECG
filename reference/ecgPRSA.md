# Phase-Rectified Signal Averaging (PRSA): DC and AC

Computes the deceleration capacity (DC) and acceleration capacity (AC)
of the RR series by phase-rectified signal averaging (Bauer et al.
2006). Anchors are points where the RR level increases (DC) or decreases
(AC); segments centered on the anchors are averaged into a PRSA curve,
and the capacity is a Haar-wavelet coefficient of that curve at scale
`s`.

## Usage

``` r
ecgPRSA(rr, T = 1L, s = 2L, L = 15L)
```

## Arguments

- rr:

  RR intervals (ms): a numeric vector or a data.frame with an `rr_ms`
  column.

- T:

  Averaging window (beats) for anchor selection (default 1: an anchor is
  a single beat larger/smaller than the previous).

- s:

  Wavelet scale for the DC/AC coefficient (default 2).

- L:

  Half-length of the averaged PRSA segments (default 15; must be \\\ge
  s\\).

## Value

A list with `dc` (deceleration capacity, ms; positive with dominant
decelerations), `ac` (acceleration capacity, ms; negative), the PRSA
curves `prsa_dc`/`prsa_ac`, the anchor counts `n_dc`/`n_ac`, and the
parameters `T`, `s`, `L`.

## References

Bauer, A. et al. (2006). Deceleration capacity of heart rate as a
predictor of mortality after myocardial infarction: cohort study.
*Lancet*, 367(9523), 1674-1681.

## See also

[`ecgHRT`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRT.md),
[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)

## Examples

``` r
set.seed(1)
rr <- 800 + cumsum(rnorm(600, sd = 3))       # RR random walk
prsa <- ecgPRSA(rr)
c(DC = prsa$dc, AC = prsa$ac)
#>        DC        AC 
#>  1.245939 -1.062157 
```
