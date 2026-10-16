function [img, info] = rda_focus(raw, t, p)
%RDA_FOCUS  Range-Doppler algorithm: range compression, azimuth FFT, RCMC, azimuth compression.
%   Zero-squint stripmap. RCMC uses a 16-tap Kaiser-windowed sinc interpolator in the range-Doppler
%   domain (an 8-tap plain sinc raised the range PSLR by about 1.5 dB; see docs/VERIFICATION.md).
c = p.c; [Na, Nr] = size(raw);
% ---- 1. range compression (matched filter in the frequency domain) ----
fr  = ((0:Nr-1) - Nr/2) * p.fs / Nr;
if isfield(p, 'range_filter') && strcmp(p.range_filter, 'posp')
    % analytic filter from the principle of stationary phase: ignores the Fresnel ripple
    Hr = exp(1j * pi * fr.^2 / p.Kr) .* (abs(fr) <= p.Br / 2);
else
    % replica filter: conjugate spectrum of the transmitted chirp, magnitude normalised,
    % so the Fresnel ripple of the real chirp spectrum is matched
    tr = ((0:Nr-1) - Nr/2) / p.fs;
    rep = (abs(tr) <= p.Tp / 2) .* exp(1j * pi * p.Kr * tr.^2);
    Rf = fftshift(fft(ifftshift(rep)));
    Hr = conj(Rf) ./ max(abs(Rf), eps) .* (abs(fr) <= p.Br / 2);
end
Hr  = Hr .* hamming_(fr / p.Br, p.win_alpha);                  % spectral weighting (1 = none)
S   = fftshift(fft(raw, [], 2), 2);
src = ifft(ifftshift(S .* repmat(Hr, Na, 1), 2), [], 2);
% the compressed pulse for a target at delay t0 peaks at fast time t0
% ---- 2. azimuth FFT ----
feta = ((0:Na-1)' - Na/2) * p.PRF / Na;
SRD  = fftshift(fft(src, [], 1), 1);
% ---- 3. RCMC ----
R  = t * c / 2;                                     % range of each range cell (m)
D  = sqrt(1 - (p.lambda * feta / (2 * p.V)).^2);    % migration factor
dR = R(:)' .* (1 ./ D - 1);                         % Na x Nr range shift (m)
dn = dR / (c / (2 * p.fs));                         % shift in range cells
if isfield(p, 'rcmc') && ~p.rcmc, dn = zeros(size(dn)); end
out = zeros(Na, Nr);
taps = -7:8; nt = numel(taps); beta = 2.5;
for a = 1:Na
    row = SRD(a, :);
    pos = (1:Nr) + dn(a, :);
    acc = zeros(1, Nr);
    base = floor(pos);
    frac = pos - base;
    for k = taps
        idx = base + k;
        x = frac - k;
        w = sinc_(x) .* kaiser_(x, nt / 2, beta);
        valid = idx >= 1 & idx <= Nr;
        acc(valid) = acc(valid) + w(valid) .* row(idx(valid));
    end
    out(a, :) = acc;
end
% ---- 4. azimuth compression ----
% Exact zero-squint filter. Only the Doppler-dependent part (D - 1) is applied: the constant
% term exp(j 4 pi fc R / c) would put a range-varying phase ramp across each target's samples
% and corrupt any band-limited interpolation of the image (for example IRF upsampling).
Ha = exp(1j * 4 * pi * R(:)' .* (D - 1) * p.fc / c);
Ha = Ha .* (abs(feta) <= p.Ba / 2);                 % process the Doppler bandwidth
Ha = Ha .* hamming_(feta / p.Ba, p.win_alpha);
img = ifft(ifftshift(out .* Ha, 1), [], 1);
info.R = R; info.feta = feta; info.max_rcm_cells = max(abs(dn(:)));
end

function y = sinc_(x)
y = ones(size(x)); nz = x ~= 0;
y(nz) = sin(pi * x(nz)) ./ (pi * x(nz));
end

function w = hamming_(u, alpha)
%HAMMING_  Generalised Hamming window over normalised frequency u in [-0.5, 0.5].
w = alpha + (1 - alpha) * cos(2 * pi * u);
w(abs(u) > 0.5) = 0;
end

function w = kaiser_(x, half, beta)
%KAISER_  Kaiser window evaluated at offsets x for a kernel of half-length 'half' (in samples).
u = x / half; w = zeros(size(x)); in = abs(u) <= 1;
w(in) = besseli(0, beta * sqrt(1 - u(in).^2)) / besseli(0, beta);
end
