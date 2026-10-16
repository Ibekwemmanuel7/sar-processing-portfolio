function p = sar_params()
%SAR_PARAMS  System parameters for a C-band spaceborne stripmap SAR (Sentinel-1-like, simplified).
%   Straight-line flight, zero squint, flat geometry with an effective velocity.
p.c      = 299792458;          % speed of light (m/s)
p.fc     = 5.405e9;            % carrier frequency (Hz)
p.lambda = p.c / p.fc;         % wavelength (m)
p.V      = 7000;               % effective platform velocity (m/s)
p.R0     = 850e3;              % slant range to scene centre (m)
p.La     = 10;                 % antenna length in azimuth (m)
p.Tp     = 10e-6;              % chirp duration (s)
p.Br     = 30e6;               % chirp bandwidth (Hz)
p.Kr     = p.Br / p.Tp;        % chirp rate (Hz/s)
p.fs     = 36e6;               % range sampling rate (Hz), 1.2 x bandwidth
p.Ba     = 0.886 * 2 * p.V / p.La;   % azimuth Doppler bandwidth (Hz)
p.PRF    = 1.25 * p.Ba;        % pulse repetition frequency (Hz)
p.Ta     = 0.886 * p.lambda * p.R0 / (p.La * p.V);   % synthetic aperture time (s)
p.Nr     = 1024;               % range samples per echo window
p.Na     = 2048;               % pulses
p.win_alpha = 0.75;           % generalised Hamming weighting (1.0 = no weighting)
p.snr_in_db = 0;               % raw-data SNR per sample for each target (dB)
% point targets: [azimuth offset (m), slant-range offset (m), amplitude]
p.targets = [   0      0   1.0;
              -150   -200   0.8;
               200    250   0.6;
               100   -350   0.5];
% example requirement (weighted processing, alpha = 0.75) and verification tolerance
p.spec.irw_rg_m  = 5.5;    % slant-range resolution, -3 dB width (m)
p.spec.irw_az_m  = 6.5;    % azimuth resolution, -3 dB width (m)
p.spec.pslr_db   = -20.0;  % peak sidelobe ratio (dB), both directions
p.spec.islr_db   = -15.0;  % integrated sidelobe ratio (dB), both directions
p.spec.loc_err_m = 1.0;    % peak location error (m)
p.spec.model_tol = 0.03;   % measured IRW must agree with the analytic model within 3 percent
p.spec.nesz_db   = -22.0;  % noise-equivalent sigma zero across the swath (dB)
end
