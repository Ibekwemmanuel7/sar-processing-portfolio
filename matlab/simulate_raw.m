function [raw, t, eta, p] = simulate_raw(p)
%SIMULATE_RAW  Simulate stripmap SAR raw echoes (baseband) for point targets.
%   Exact hyperbolic range history, LFM chirp, two-way sinc^2 azimuth antenna pattern,
%   complex white Gaussian noise at p.snr_in_db per sample.
c = p.c;
Rmin = p.R0 - 0.5 * (p.Nr / p.fs) * c / 2 + p.Tp * c / 4;   % window start, centres the scene
t   = 2 * Rmin / c + (0:p.Nr-1) / p.fs;            % fast time (s)
eta = ((0:p.Na-1) - p.Na/2) / p.PRF;               % slow time (s)
[T, ETA] = meshgrid(t, eta);                        % Na x Nr
raw = zeros(p.Na, p.Nr);
theta_bw = 0.886 * p.lambda / p.La;                 % azimuth beamwidth (rad)
for k = 1:size(p.targets, 1)
    xa = p.targets(k, 1); dr = p.targets(k, 2); A = p.targets(k, 3);
    R0k = p.R0 + dr;
    R = sqrt(R0k^2 + (p.V * ETA - xa).^2);         % range history
    tau = T - 2 * R / c;
    theta = atan((p.V * ETA - xa) / R0k);
    wa = sinc_(0.886 * theta / theta_bw).^2;        % two-way antenna pattern
    wr = abs(tau) <= p.Tp / 2;                      % chirp envelope
    raw = raw + A * wr .* wa .* exp(-1j * 4 * pi * p.fc * R / c) .* exp(1j * pi * p.Kr * tau.^2);
end
sig_pow = 1;                                        % unit-amplitude reference target
noise_std = sqrt(sig_pow / 10^(p.snr_in_db / 10) / 2);
randn('seed', 7);
raw = raw + noise_std * (randn(size(raw)) + 1j * randn(size(raw)));
end

function y = sinc_(x)
y = ones(size(x)); nz = x ~= 0;
y(nz) = sin(pi * x(nz)) ./ (pi * x(nz));
end
