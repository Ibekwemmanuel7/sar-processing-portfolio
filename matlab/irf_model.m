function mdl = irf_model(p)
%IRF_MODEL  Analytic (model) impulse response from the processed spectra, for verification.
%   Range: rectangular band Br with generalised Hamming weighting.
%   Azimuth: band Ba weighted by the two-way antenna pattern mapped to Doppler and by Hamming.
N = 2^16;
% range
f = ((0:N-1) - N/2) / N * 8 * p.Br; df = f(2) - f(1);
W = (abs(f) <= p.Br/2) .* (p.win_alpha + (1 - p.win_alpha) * cos(2*pi*f/p.Br));
[mdl.irw_rg_m, mdl.pslr_rg_db, mdl.islr_rg_db] = irf_from_spectrum(W, df, p.c / 2);
% azimuth
f = ((0:N-1) - N/2) / N * 8 * p.Ba; df = f(2) - f(1);
theta = p.lambda * f / (2 * p.V); theta_bw = 0.886 * p.lambda / p.La;
x = 0.886 * theta / theta_bw; ant = ones(size(x)); nz = x ~= 0; ant(nz) = (sin(pi*x(nz)) ./ (pi*x(nz))).^2;
W = ant .* (abs(f) <= p.Ba/2) .* (p.win_alpha + (1 - p.win_alpha) * cos(2*pi*f/p.Ba));
[mdl.irw_az_m, mdl.pslr_az_db, mdl.islr_az_db] = irf_from_spectrum(W, df, p.V);
end

function [irw, pslr, islr] = irf_from_spectrum(W, df, scale)
N = numel(W);
h = abs(fftshift(ifft(ifftshift(W)))).^2; h = h / max(h);
[~, k0] = max(h);
dt = 1 / (N * df);                      % time (or slow-time) spacing
[irw, pslr, islr] = cut_metrics(h, k0, dt * scale);
end
