# Simulate PPG with known fiducials

Generates DC-offset PPG pulses from systolic and diastolic Gaussian
waves plus a sharp foot dip, with optional respiratory amplitude
modulation.

## Usage

``` r
make_ppg(
  n_time = 15000,
  n_channels = 1,
  sr = 125,
  heart_rate = 72,
  dc = 2,
  sys_amp = 1,
  dia_amp = 0.4,
  crest_ms = 150,
  dia_offset_ms = 350,
  noise_sd = 0.02,
  foot_dip = 0.1,
  resp_freq = NULL,
  resp_depth = 0.3
)
```

## Arguments

- n_time:

  Number of samples.

- n_channels:

  Number of PPG channels.

- sr:

  Sampling rate in Hz.

- heart_rate:

  Pulse rate in beats per minute.

- dc:

  Non-pulsatile baseline level.

- sys_amp, dia_amp:

  Systolic and diastolic wave amplitudes.

- crest_ms:

  Foot-to-systolic interval in milliseconds.

- dia_offset_ms:

  Foot-to-diastolic interval in milliseconds.

- noise_sd:

  Gaussian noise standard deviation.

- foot_dip:

  Depth of the local foot minimum.

- resp_freq:

  Optional respiratory amplitude-modulation frequency in Hz.

- resp_depth:

  Fractional respiratory modulation depth.

## Value

A list with `pe`, a `PhysioExperiment`, and `truth`, a fiducial data
frame (`beat`, `foot`, `systolic`, `diastolic`).

## Examples

``` r
set.seed(1)
simulation <- make_ppg(n_time = 3000)
head(simulation$truth)
#>   beat foot systolic diastolic
#> 1    1  104      123       148
#> 2    2  208      227       252
#> 3    3  312      331       356
#> 4    4  416      435       460
#> 5    5  520      539       564
#> 6    6  624      643       668
```
