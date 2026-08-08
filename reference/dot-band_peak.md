# Peak frequency within a band

Peak frequency within a band

## Usage

``` r
.band_peak(freqs, psd, band)
```

## Arguments

- freqs:

  Frequency vector (Hz).

- psd:

  Power spectral density vector.

- band:

  Numeric length-2 band (Hz).

## Value

Frequency (Hz) of maximum PSD within the band, or NA if empty.
