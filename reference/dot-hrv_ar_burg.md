# Burg autoregressive PSD for HRV

Burg autoregressive PSD for HRV

## Usage

``` r
.hrv_ar_burg(t_sec, rr_ms, order = NULL, detrend = FALSE, lambda = 500)
```

## Arguments

- t_sec:

  Time stamps in seconds.

- rr_ms:

  RR intervals in milliseconds.

- order:

  AR model order, or NULL to select by AIC.

- detrend:

  Logical; apply smoothness-priors detrending after resampling.

- lambda:

  Smoothing parameter for smoothness-priors detrending.

## Value

List with freqs and psd vectors.
