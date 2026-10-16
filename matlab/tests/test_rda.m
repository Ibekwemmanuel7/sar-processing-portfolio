%TEST_RDA  Regression tests for the simulator and range-Doppler processor (MATLAB or Octave).
%   Run from the matlab folder:  run('tests/test_rda.m')
addpath(fileparts(fileparts(mfilename('fullpath'))));
p = sar_params(); p.targets = [0 0 1]; p.snr_in_db = 40;
[raw, t, eta, p] = simulate_raw(p);
% 1. raw data is not focused: energy is spread over many cells
assert(max(abs(raw(:))) < 2, 'raw data should not contain a focused peak');
% 2. weighted processing meets the requirement on an isolated target
p.range_filter = 'replica';
[img, info] = rda_focus(raw, t, p);
m = measure_irf(img, p, t, eta, 0, 0); mdl = irf_model(p); s = p.spec;
assert(m.irw_rg_m <= s.irw_rg_m && m.irw_az_m <= s.irw_az_m, 'resolution requirement');
assert(max(m.pslr_rg_db, m.pslr_az_db) <= s.pslr_db, 'PSLR requirement');
assert(max(m.islr_rg_db, m.islr_az_db) <= s.islr_db, 'ISLR requirement');
assert(abs(m.irw_rg_m / mdl.irw_rg_m - 1) <= s.model_tol, 'range IRW agrees with the model');
assert(abs(m.irw_az_m / mdl.irw_az_m - 1) <= s.model_tol, 'azimuth IRW agrees with the model');
assert(max(abs([m.loc_err_rg_m m.loc_err_az_m])) <= s.loc_err_m, 'peak location');
% 3. unweighted processing reproduces the sinc sidelobe level within 1 dB
p.win_alpha = 1.0; img1 = rda_focus(raw, t, p); m1 = measure_irf(img1, p, t, eta, 0, 0);
assert(abs(m1.pslr_rg_db - (-13.26)) < 1.0, 'unweighted range PSLR near -13.3 dB');
% 4. a plain 8-tap sinc RCMC interpolator was found to raise range PSLR; the Kaiser kernel must hold it
%    within 1 dB of the analytic model
assert(m.pslr_rg_db <= mdl.pslr_rg_db + 1.0, 'RCMC interpolation keeps range PSLR within 1 dB of model');
fprintf('test_rda: all assertions passed\n');
