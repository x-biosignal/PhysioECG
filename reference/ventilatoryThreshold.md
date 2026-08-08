# Detect ventilatory thresholds from incremental exercise data

VT1 is estimated by V-slope breakpoint regression of VCO2 on VO2. VT2 is
estimated above VT1 from the breakpoint of VE on VCO2.

## Usage

``` r
ventilatoryThreshold(
  vo2,
  vco2 = NULL,
  ve = NULL,
  time = NULL,
  average_n = 1L,
  min_segment = 3L
)
```

## Arguments

- vo2:

  Numeric VO2 vector, or a data frame containing `vo2`, `vco2`, and
  `ve`.

- vco2, ve:

  Numeric vectors matching `vo2` when it is not a data frame.

- time:

  Optional time vector. For data-frame input, an existing `time` column
  is used unless this argument is supplied.

- average_n:

  Positive integer number of ordered samples per block mean.

- min_segment:

  Minimum number of observations fitted on each side of a breakpoint.

## Value

A `ventilatory_threshold` object with VT1, VT2, and the data actually
fitted.

## References

Beaver WL, Wasserman K, Whipp BJ (1986). A new method for detecting
anaerobic threshold by gas exchange. *Journal of Applied Physiology*,
60:2020-2027.
[doi:10.1152/jappl.1986.60.6.2020](https://doi.org/10.1152/jappl.1986.60.6.2020)

## Examples

``` r
vo2 <- seq(1, 3, by = 0.05)
vco2 <- ifelse(vo2 <= 2, 0.9 * vo2, 1.8 + 1.3 * (vo2 - 2))
ventilatoryThreshold(vo2, vco2, ve = 25 * vco2)
#> <ventilatory_threshold> method=v-slope
#>   VT1: VO2 2, VCO2 1.8, VE 45 (slopes 0.9 -> 1.3)
#>   VT2: not detected
```
