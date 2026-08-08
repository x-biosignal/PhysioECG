# Welch's method for HRV PSD estimation

Welch's method for HRV PSD estimation

## Usage

``` r
.hrv_welch(t_sec, rr_ms, detrend = FALSE, lambda = 500)
```

## Arguments

- t_sec:

  Time stamps in seconds.

- rr_ms:

  RR intervals in milliseconds.

- detrend:

  Logical; apply smoothness-priors detrending after resampling.

- lambda:

  Smoothing parameter for smoothness-priors detrending.

## Value

List with freqs and psd vectors.
