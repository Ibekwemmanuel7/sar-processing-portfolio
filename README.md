# SAR processing portfolio: simulator, range-Doppler processor and verification

Emmanuel Chidi Ibekwe · M.Sc. Atmospheric Science, Texas A&M University · Houston, TX

This repository shows end-to-end SAR processing work in **MATLAB** (runs in MATLAB or GNU Octave) with a
matching **Python** implementation used as an independent cross-check:

| Job requirement | Where it is in this repository |
|---|---|
| Develop algorithms to process SAR data | `matlab/rda_focus.m`: range-Doppler algorithm (replica range compression, azimuth FFT, RCMC with a Kaiser-windowed sinc interpolator, exact zero-squint azimuth compression, spectral weighting). `detect_dark.m`: adaptive-threshold oil-slick detector. |
| Model and verify system performance to meet specifications | `irf_model.m` (analytic IRF from the processed spectra), `measure_irf.m` (IRW, PSLR, ISLR, location error), `nesz_model.m` (radar-equation NESZ), `run_all.m` (requirement checks). See `docs/VERIFICATION.md`. |
| Develop and maintain simulators for testing and verification | `simulate_raw.m` (point-target raw-data simulator: LFM chirp, hyperbolic range history, two-way antenna pattern, noise). `slick_scene.m` (K-distributed sea clutter with an oil slick and a low-wind look-alike, with truth masks). Regression tests in `matlab/tests` and `python/tests`. |
| Collaborate with ground and software teams | Same algorithms in two languages with agreeing results, scripted runs, CSV outputs, tests that run headless (CI-ready). |

![raw and focused](docs/figures/fig1_raw_vs_focused.png)

## Results (C-band stripmap, Sentinel-1-like parameters)

Operational case: Hamming weighting (alpha = 0.75), replica range filter. Four point targets.

| Metric | Requirement | Analytic model | Measured (4 targets, MATLAB/Octave) | Result |
|---|---|---|---|---|
| Slant-range IRW | <= 5.5 m | 4.99 m | 5.03 to 5.04 m | PASS |
| Azimuth IRW | <= 6.5 m | 6.21 m | 6.24 to 6.27 m | PASS |
| Measured IRW vs model | within 3 % | | +0.8 to 1.0 % (range), +0.5 to 1.0 % (azimuth) | PASS |
| PSLR (range / azimuth) | <= -20 dB | -21.4 / -26.9 dB | -20.5 to -20.9 / -26.3 to -27.5 dB | PASS |
| ISLR (range / azimuth) | <= -15 dB | -16.0 / -22.0 dB | -16.3 to -16.4 / -22.5 to -22.9 dB | PASS |
| Peak location error | <= 1 m | 0 | <= 0.14 m | PASS |
| NESZ, 20 to 46 deg incidence | <= -22 dB | -28.8 to -22.2 dB | (model) | PASS, 0.2 dB margin at far range |

The unweighted case reproduces the sinc response (range PSLR -13.1 dB against -13.3 dB theory).
The Python implementation gives the same IRW to within 0.01 m, and PSLR to within about 1 dB (different noise realisation).

![IRF cuts](docs/figures/fig2_irf_cuts.png)

### A defect found by the verification
The first processor used an 8-tap plain sinc interpolator for range cell migration correction. The
verification run showed range PSLR of -19.7 to -20.2 dB, so three of four targets failed the -20 dB
requirement, while the analytic model predicted -21.4 dB. Disabling RCMC on an isolated target
restored -21.7 dB, which isolated the interpolator. A 16-tap Kaiser-windowed sinc kernel brings the
range PSLR to -20.5 to -20.9 dB. A regression test now holds the range PSLR within 1 dB of the model.
A second defect: applying the full azimuth filter phase exp(j4πfcR/c) per range cell put a
range-varying phase ramp on each target and distorted the upsampled response (range IRW measured
8 % below theory). The filter now applies only the Doppler-dependent term.

## Oil-slick detection against simulated truth

![slick](docs/figures/fig4_slick_simulation.png)

Probability of detection rises from 0 at 2 to 3 dB contrast to 0.80 at 6 dB and 0.91 at 10 dB, with
open-sea false-alarm rates of 1e-5 to 2e-4. The broad low-wind look-alike is not flagged, because the
adaptive background window follows it; a look-alike smaller than the window would be. This is why
dark-spot detection is a first stage, followed by shape and context checks against a pre-event image.

## Applied work on real Sentinel-1 data

![Sentinel-1 case](docs/figures/fig5_sentinel1_case_study.png)

Sentinel-1 GRD time series (Google Earth Engine) over an FPSO that exploded off Delta State, Nigeria,
on 2 February 2022: the vessel is a bright point target before the event, and dark slicks trail from
its position 6 and 18 days after. Fire detections (NASA FIRMS) confirm the event date. Other case
studies: the MV Wakashio spill (Mauritius, 2020) and onshore pipeline-spill monitoring in the Niger
Delta (radar change detection, Sentinel-2 indices, VIIRS fire data).

## How to run

MATLAB or Octave, from `matlab/`:

    run_all              % simulate, focus, measure, verify; writes results/
    run_slick            % slick detector sweep over contrast and looks
    run('tests/test_rda.m')

Python 3.10+, from `python/`:

    pip install numpy scipy matplotlib pytest
    python run_all.py    # figures in docs/figures
    python -m pytest -q  # 5 tests

Runtime: about 20 s (Octave) or 10 s (Python) on a laptop. Raw data 2048 x 1024 complex samples.

## Simplifications (stated so the results are read correctly)
Straight-line geometry with an effective velocity, zero squint, no Earth curvature or Doppler
centroid estimation, no secondary range compression (not needed at this squint and bandwidth),
point targets only for the IRF tests, and an image-domain (not raw-domain) clutter simulator for the
slick tests. Natural extensions: squinted geometry and secondary range compression, chirp scaling or
omega-k, Doppler centroid and autofocus, and ambiguity (AASR, RASR) analysis.

## References
Cumming, I. G. and Wong, F. H. (2005). *Digital Processing of Synthetic Aperture Radar Data.* Artech House.
Curlander, J. C. and McDonough, R. N. (1991). *Synthetic Aperture Radar: Systems and Signal Processing.* Wiley.

MIT licence. Simulated data only, except the Sentinel-1 figure: contains modified Copernicus Sentinel data (2022).
