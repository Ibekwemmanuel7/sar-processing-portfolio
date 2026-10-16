# Verification record

Run: `matlab/run_all.m` in GNU Octave (same code runs in MATLAB). Full output: `octave_verification.txt`.

## Method
1. Simulate raw echoes for four point targets (input SNR 0 dB per sample).
2. Focus with the range-Doppler algorithm, three cases: unweighted with replica range filter; Hamming
   (alpha 0.75) with the analytic POSP filter; Hamming with the replica filter (operational).
3. For each target, upsample a 32 x 32 patch by 16 and measure -3 dB width (IRW), peak sidelobe ratio
   (PSLR), integrated sidelobe ratio (ISLR, sidelobes outside the first nulls) and peak location error.
4. Compute the analytic IRF from the processed spectra (range: rectangle x Hamming; azimuth: two-way
   antenna pattern mapped to Doppler x rectangle x Hamming) and compare.
5. Check each metric against the requirement table in `sar_params.m`, and the measured IRW against the
   model within 3 percent. Compute NESZ across the swath from the radar equation.

## Requirement table (example, weighted processing)

| Parameter | Requirement |
|---|---|
| Slant-range IRW | <= 5.5 m |
| Azimuth IRW | <= 6.5 m |
| PSLR, both directions | <= -20 dB |
| ISLR, both directions | <= -15 dB |
| Peak location error | <= 1 m |
| IRW agreement with the model | within 3 % |
| NESZ, 20 to 46 deg incidence | <= -22 dB |

## Findings
- **RCMC interpolator (fixed).** 8-tap plain sinc: range PSLR -19.7 to -20.2 dB (3 of 4 targets fail).
  Isolated by disabling RCMC on a single high-SNR target (-21.7 dB). 16-tap Kaiser-windowed sinc
  (beta 2.5): -20.5 to -20.9 dB. Remaining gap to the model (0.5 to 0.9 dB) comes from interpolating a
  signal sampled at only 1.2 times its bandwidth.
- **Azimuth filter phase (fixed).** Applying exp(j 4 pi fc R / c) per range cell put a fast phase ramp
  across each target's range samples, and the upsampled range IRW read 4.07 m against 4.43 m theory.
  Applying only the Doppler-dependent term (D - 1) gives 4.46 m (unweighted), within 2.3 % of the model.
- **Range filter.** Replica and POSP filters give the same result here (time-bandwidth product 300), so
  the PSLR gap was not caused by the Fresnel ripple of the chirp spectrum.
- **NESZ margin.** The far-range NESZ (-22.2 dB at 46 deg) has only 0.2 dB margin. Options: more peak
  power, a longer pulse, or limiting the swath to 44 deg.
