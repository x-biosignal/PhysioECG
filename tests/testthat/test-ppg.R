library(testthat)


test_that("ppgDetectPulses recovers systolic peaks with high Se and PPV", {
  set.seed(1)
  simulation <- make_ppg(
    n_time = 15000, sr = 125, heart_rate = 72
  )
  peaks <- ppgDetectPulses(simulation$pe)

  expect_identical(
    names(peaks),
    c("channel", "sample", "time_sec", "amplitude")
  )
  tolerance <- as.integer(round(0.15 * 125))
  true_positive <- PhysioECG:::.ecg_match_beats(
    sort(peaks$sample),
    sort(simulation$truth$systolic),
    tolerance
  )
  sensitivity <- true_positive / nrow(simulation$truth)
  ppv <- true_positive / nrow(peaks)
  expect_gt(sensitivity, 0.95)
  expect_gt(ppv, 0.95)
  expect_equal(
    peaks$amplitude,
    SummarizedExperiment::assay(simulation$pe)[peaks$sample, 1]
  )
})


test_that("ppgDetectPulses compensates FIR group delay", {
  set.seed(11)
  simulation <- make_ppg(n_time = 5000, noise_sd = 0)
  peaks <- ppgDetectPulses(simulation$pe)
  tolerance <- as.integer(round(0.08 * 125))
  nearest_error <- vapply(
    simulation$truth$systolic,
    function(truth) min(abs(peaks$sample - truth)),
    integer(1)
  )

  expect_lt(stats::median(nearest_error), 2)
  expect_true(mean(nearest_error <= tolerance) > 0.95)
})


test_that("ppgDetectPulses returns typed empty output for silence", {
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = matrix(2, 1000, 1)),
    colData = S4Vectors::DataFrame(label = "PPG1", type = "PPG"),
    samplingRate = 125
  )
  peaks <- ppgDetectPulses(experiment)
  expect_equal(nrow(peaks), 0L)
  expect_type(peaks$channel, "integer")
  expect_type(peaks$sample, "integer")
})


test_that("ppgDetectPulses handles short and high-offset records", {
  short <- PhysioCore::PhysioExperiment(
    assays = list(raw = matrix(seq_len(20), ncol = 1)),
    colData = S4Vectors::DataFrame(label = "PPG1", type = "PPG"),
    samplingRate = 125
  )
  expect_equal(nrow(ppgDetectPulses(short)), 0L)

  sr <- 100
  time <- seq(0, 10 - 1 / sr, by = 1 / sr)
  signal <- 1e9 + sin(2 * pi * 1.2 * time)
  offset <- PhysioCore::PhysioExperiment(
    assays = list(raw = matrix(signal, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "PPG1", type = "PPG"),
    samplingRate = sr
  )
  expect_gt(nrow(ppgDetectPulses(offset)), 5L)
  expect_true(is.finite(ppgQuality(offset)$s_sqi))
})


test_that("PRV time and frequency equal HRV under constant PTT", {
  set.seed(2)
  ecg <- make_ecg(n_time = 60000, sr = 500, heart_rate = 72)
  r_peaks <- ecgDetectRpeaks(ecg)
  rr <- ecgRRintervals(ecg, r_peaks)
  hrv_time <- ecgHRVtime(rr, rhythm_check = FALSE)
  hrv_frequency <- ecgHRVfreq(rr, rhythm_check = FALSE)

  ppg_peaks <- r_peaks
  ppg_peaks$sample <- r_peaks$sample +
    as.integer(round(0.250 * 500))
  ppg_peaks$time_sec <- r_peaks$time_sec
  prv <- pulseRateVariability(ecg, peaks = ppg_peaks, freq = TRUE)

  expect_equal(prv$time$sdpp, hrv_time$sdnn, tolerance = 1e-6)
  expect_equal(prv$time$rmssd, hrv_time$rmssd, tolerance = 1e-6)
  expect_equal(prv$time$mean_pp, hrv_time$mean_rr, tolerance = 1e-6)
  expect_equal(prv$freq$lf, hrv_frequency$lf, tolerance = 1e-6)
  expect_equal(prv$freq$hf, hrv_frequency$hf, tolerance = 1e-6)
  expect_false(any(c("rhythm", "hrv_valid") %in% names(prv$time)))
})


test_that("PRV variance follows the PTT-jitter identity", {
  set.seed(3)
  n <- 4000
  jitter_sd <- 12
  rr_ms <- 800 + 40 * sin(2 * pi * 0.1 * seq_len(n)) +
    rnorm(n, 0, 20)
  transit <- rnorm(n + 1L, 0, jitter_sd)
  ppi_ms <- rr_ms + diff(transit)
  rmssd <- function(values) sqrt(mean(diff(values)^2))

  expect_lt(
    abs((stats::var(ppi_ms) - stats::var(rr_ms)) -
          2 * jitter_sd^2),
    0.15 * (2 * jitter_sd^2)
  )
  expect_lt(
    abs((rmssd(ppi_ms)^2 - rmssd(rr_ms)^2) -
          6 * jitter_sd^2),
    0.25 * (6 * jitter_sd^2)
  )
})


test_that("pulseRateVariability handles empty and optional domains", {
  simulation <- make_ppg(n_time = 3000, noise_sd = 0)
  empty <- data.frame(
    channel = integer(),
    sample = integer(),
    time_sec = numeric(),
    amplitude = numeric()
  )
  result <- pulseRateVariability(
    simulation$pe, peaks = empty, freq = FALSE, nonlinear = FALSE
  )
  expect_equal(nrow(result$ppi), 0L)
  expect_equal(nrow(result$time), 0L)
  expect_null(result$freq)
  expect_null(result$nonlinear)
})


test_that("perfusionIndex equals closed-form AC over DC", {
  sr <- 100
  time <- seq(0, 20 - 1 / sr, by = 1 / sr)
  ppg <- 2 + 0.1 * sin(2 * pi * 1.2 * time)
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = matrix(ppg, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "PPG1", type = "PPG"),
    samplingRate = sr
  )
  peaks <- ppgDetectPulses(experiment)
  result <- perfusionIndex(experiment, peaks = peaks)

  expect_lt(abs(result$dc - 2), 0.01)
  expect_lt(abs(result$ac - 0.2), 0.02)
  expect_lt(abs(result$pi_pct - 10), 0.5)
  expect_equal(result$n_pulses, nrow(peaks))
})


test_that("perfusionIndex guards absent pulses and zero DC", {
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = cbind(flat = rep(2, 100), zero = rep(0, 100))),
    colData = S4Vectors::DataFrame(
      label = c("flat", "zero"), type = c("PPG", "PPG")
    ),
    samplingRate = 20
  )
  empty <- data.frame(channel = integer(), sample = integer())
  result <- perfusionIndex(experiment, peaks = empty)
  expect_equal(result$n_pulses, c(0L, 0L))
  expect_true(all(is.na(result$ac)))
  expect_true(all(is.na(result$pi_pct)))
})


test_that("pulseWaveFeatures recovers rise time, AI, RI, PTT, and PWV", {
  set.seed(4)
  simulation <- make_ppg(
    n_time = 8000,
    sr = 125,
    heart_rate = 60,
    sys_amp = 1,
    dia_amp = 0.4,
    crest_ms = 160
  )
  peaks <- data.frame(
    channel = 1L,
    sample = simulation$truth$systolic,
    time_sec = (simulation$truth$systolic - 1) / 125,
    amplitude = NA_real_
  )
  ecg_peaks <- data.frame(
    sample = simulation$truth$foot -
      as.integer(round(0.200 * 125))
  )
  features <- pulseWaveFeatures(
    simulation$pe,
    peaks = peaks,
    ecg_peaks = ecg_peaks,
    distance_m = 0.5
  )
  middle <- features[3:(nrow(features) - 1L), ]

  expect_lt(abs(stats::median(middle$rise_time_ms) - 160), 8)
  expect_lt(
    abs(stats::median(middle$augmentation_index, na.rm = TRUE) - 54.5),
    5
  )
  expect_lt(
    abs(stats::median(middle$reflection_index, na.rm = TRUE) - 45.5),
    5
  )
  expect_equal(
    middle$augmentation_index + middle$reflection_index,
    rep(100, nrow(middle)),
    tolerance = 1e-10
  )
  expect_lt(abs(stats::median(middle$ptt_ms, na.rm = TRUE) - 200), 8)
  expect_lt(abs(stats::median(middle$pwv_m_s, na.rm = TRUE) - 2.5), 0.3)
})


test_that("pulseWaveFeatures uses channel-matched preceding ECG peaks", {
  simulation <- make_ppg(
    n_time = 3000, n_channels = 2, noise_sd = 0
  )
  peaks <- do.call(rbind, lapply(1:2, function(channel) {
    data.frame(
      channel = channel,
      sample = simulation$truth$systolic,
      time_sec = (simulation$truth$systolic - 1) / 125
    )
  }))
  ecg_peaks <- rbind(
    data.frame(channel = 1L, sample = simulation$truth$foot - 25L),
    data.frame(channel = 2L, sample = simulation$truth$foot - 30L)
  )
  features <- pulseWaveFeatures(
    simulation$pe, peaks, ecg_peaks, distance_m = 0.5
  )

  expect_equal(
    unique(features$ptt_ms[features$channel == 1L]), 200
  )
  expect_equal(
    unique(features$ptt_ms[features$channel == 2L]), 240
  )
})


test_that("ppgQuality matches population moment calculations", {
  values <- c(0, 0, 0, 4)
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = matrix(values, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "PPG1", type = "PPG"),
    samplingRate = 1
  )
  quality <- ppgQuality(experiment)

  expect_equal(quality$s_sqi, 24 / (4 * 3^1.5), tolerance = 1e-6)
  expect_equal(quality$k_sqi, 84 / 36, tolerance = 1e-6)
  expect_true(quality$accept)

  windowed <- ppgQuality(
    experiment, window_sec = 2, skew_threshold = 0.1
  )
  expect_equal(nrow(windowed), 2L)
  expect_equal(windowed$start_sec, c(0, 2))
  expect_false(any(windowed$accept))
})


test_that("ppgQuality marks constant signals unusable", {
  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = matrix(2, 100, 1)),
    colData = S4Vectors::DataFrame(label = "PPG1", type = "PPG"),
    samplingRate = 10
  )
  quality <- ppgQuality(experiment)
  expect_true(is.na(quality$s_sqi))
  expect_true(is.na(quality$k_sqi))
  expect_false(quality$accept)
})


test_that("ppgRespiration recovers amplitude-modulation frequency", {
  set.seed(5)
  simulation <- make_ppg(
    n_time = 30000,
    sr = 125,
    heart_rate = 72,
    resp_freq = 0.25,
    resp_depth = 0.3
  )
  peaks <- ppgDetectPulses(simulation$pe)
  respiration <- ppgRespiration(
    simulation$pe, peaks = peaks, method = "am"
  )

  expect_identical(respiration$method, "am")
  expect_false(is.null(respiration$edr))
  expect_lt(
    abs(respiration$resp_rate$resp_rate_hz[1] - 0.25),
    0.05
  )
})


test_that("ppgRespiration supports BW/FM and flat-feature guards", {
  simulation <- make_ppg(
    n_time = 8000, noise_sd = 0, resp_freq = NULL
  )
  peaks <- data.frame(
    channel = 1L,
    sample = simulation$truth$systolic,
    time_sec = (simulation$truth$systolic - 1) / 125
  )
  bw <- ppgRespiration(simulation$pe, peaks, method = "bw")
  fm <- ppgRespiration(simulation$pe, peaks, method = "fm")

  expect_true(is.na(bw$resp_rate$resp_rate_hz))
  expect_true(is.na(fm$resp_rate$resp_rate_hz))
  expect_true(all(c("channel", "time_sec", "feature") %in% names(fm$beats)))
})


test_that("make_ppg exposes exact fiducial geometry", {
  set.seed(26)
  simulation <- make_ppg(
    n_time = 4000,
    n_channels = 2,
    sr = 125,
    heart_rate = 60,
    crest_ms = 160,
    dia_offset_ms = 360,
    noise_sd = 0
  )
  assay <- SummarizedExperiment::assay(simulation$pe)
  expect_equal(dim(assay), c(4000, 2))
  expect_equal(
    simulation$truth$systolic - simulation$truth$foot,
    rep(20L, nrow(simulation$truth))
  )
  expect_equal(
    simulation$truth$diastolic - simulation$truth$foot,
    rep(45L, nrow(simulation$truth))
  )
})


test_that("PPG APIs reject invalid sampling and peak metadata", {
  simulation <- make_ppg(n_time = 3000)
  expect_error(
    ppgDetectPulses(simulation$pe, high = 100),
    "Nyquist"
  )
  expect_error(
    perfusionIndex(
      simulation$pe,
      data.frame(channel = 2L, sample = 100L)
    ),
    "invalid channel"
  )
  expect_error(
    pulseWaveFeatures(
      simulation$pe,
      peaks = data.frame(channel = 1L, sample = 100L),
      distance_m = 0
    ),
    "distance_m"
  )
  expect_error(
    pulseWaveFeatures(
      simulation$pe,
      peaks = data.frame(channel = 1L, sample = 100L),
      ecg_peaks = data.frame(sample = 0L)
    ),
    "positive"
  )
  expect_error(
    make_ppg(resp_depth = 1),
    "less than 1"
  )
})
