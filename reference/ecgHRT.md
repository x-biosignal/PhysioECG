# Heart-Rate Turbulence (HRT)

Computes turbulence onset (TO) and turbulence slope (TS), the biphasic
sinus-rhythm response to a ventricular premature beat (Schmidt et al.
1999). RR intervals are indexed by beat so that `rr[i]` is the interval
between beat `i` and beat `i+1`; for a PVC at beat `p` the coupling
interval is `rr[p-1]` and the compensatory pause `rr[p]`, which are
excluded. The local tachograms of all usable PVCs are averaged before
TO/TS.

## Usage

``` r
ecgHRT(rr, pvc_index, n_post = 15L)
```

## Arguments

- rr:

  RR intervals (ms): a numeric vector or a data.frame with an `rr_ms`
  column.

- pvc_index:

  Integer beat indices of the ventricular premature beats (e.g. the
  `beat` column of
  [`ecgDetectPVC`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectPVC.md)
  where `is_pvc`).

- n_post:

  Number of post-PVC sinus intervals used for TS (default 15).

## Value

A list with `to` (turbulence onset, percent; negative is normal), `ts`
(turbulence slope, ms per RR interval; positive is normal), `n_pvc`
(PVCs averaged) and `tachogram` (the averaged post-PVC sinus RR series).

## References

Schmidt, G. et al. (1999). Heart-rate turbulence after ventricular
premature beats as a predictor of mortality after acute myocardial
infarction. *Lancet*, 353(9162), 1390-1396.

## See also

[`ecgDetectPVC`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectPVC.md),
[`ecgPRSA`](https://x-biosignal.github.io/PhysioECG/reference/ecgPRSA.md)

## Examples

``` r
# canonical post-PVC pattern: early acceleration then deceleration
rr <- c(rep(800, 5), 500, 1100, 780, 785, 790 + 5 * (0:12))
ecgHRT(rr, pvc_index = 7)
#> $to
#> [1] -2.1875
#> 
#> $ts
#> [1] 5
#> 
#> $n_pvc
#> [1] 1
#> 
#> $tachogram
#>  [1] 780 785 790 795 800 805 810 815 820 825 830 835 840 845 850
#> 
#> $to_per_pvc
#> [1] -2.1875
#> 
```
