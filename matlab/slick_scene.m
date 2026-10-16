function [I, truth, lookalike, s] = slick_scene(s)
%SLICK_SCENE  Image-domain sea-clutter simulator with an oil slick and a low-wind look-alike.
%   Intensity I = sigma0 x texture x speckle, where
%     speckle  ~ Gamma(L, 1/L)       (L-look intensity, fully developed)
%     texture  ~ Gamma(nu, 1/nu)     (K-distributed sea clutter; nu = Inf gives no texture)
%   The slick multiplies sigma0 by 10^(-contrast_db/10) inside a random elongated patch (truth).
%   The look-alike is a broad, smooth low-wind area with weaker damping (not oil).
if nargin < 1, s = struct(); end
d = struct('n', 512, 'dx', 10, 'sigma0_db', -15, 'looks', 4, 'nu', 8, 'contrast_db', 6, ...
           'lookalike_db', 4.5, 'seed', 1);
f = fieldnames(d); for k = 1:numel(f), if ~isfield(s, f{k}), s.(f{k}) = d.(f{k}); end, end
rand('seed', s.seed); randn('seed', s.seed);
n = s.n; [X, Y] = meshgrid(1:n, 1:n);
% slick: rotated ellipse with a ragged edge (smoothed noise), about 3 km x 0.6 km
th = 35 * pi / 180; cx = 0.45 * n; cy = 0.5 * n; a = 150; b = 30;
u = ((X - cx) * cos(th) + (Y - cy) * sin(th)) / a;
v = (-(X - cx) * sin(th) + (Y - cy) * cos(th)) / b;
edge = conv2(randn(n), ones(15) / 225, 'same') * 6;
truth = (u.^2 + v.^2 + 0.35 * edge) < 1;
% look-alike: broad low-wind zone in a corner, smooth Gaussian damping profile
la = exp(-(((X - 0.82 * n) / 60).^2 + ((Y - 0.2 * n) / 45).^2));
lookalike = la > 0.5;
damp_db = s.contrast_db * truth + s.lookalike_db * la;
sig = 10^(s.sigma0_db / 10) * 10.^(-damp_db / 10);
if isinf(s.nu), tex = ones(n); else, tex = gamrnd_(s.nu, n) / s.nu; end
spk = gamrnd_(s.looks, n) / s.looks;
I = sig .* tex .* spk;
end

function g = gamrnd_(k, n)
%GAMRND_  Gamma(k, 1) samples by summing exponentials for integer k (Octave/MATLAB, no toolbox).
k = round(k); g = zeros(n);
for i = 1:k, g = g - log(rand(n)); end
end
