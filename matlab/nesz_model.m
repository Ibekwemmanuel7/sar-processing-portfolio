function T = nesz_model(p)
%NESZ_MODEL  Noise-equivalent sigma zero across the swath from the radar equation.
%   NESZ = 2 (4 pi)^3 R^3 V k T0 F Br L sin(theta_i) / (Pavg G^2 lambda^3 c),  Pavg = Pt Tp PRF
%   (Curlander and McDonough, 1991). Returns [incidence deg, slant range km, NESZ dB].
k = 1.380649e-23; T0 = 290;
Pt = 4000;             % peak transmit power (W)
Wa = 1.0;              % antenna height (m)
eta_ant = 0.6;         % aperture efficiency
G  = 4 * pi * eta_ant * p.La * Wa / p.lambda^2;
F  = 10^(3.5 / 10);    % receiver noise figure 3.5 dB
L  = 10^(4.0 / 10);    % system losses 4 dB
H  = 693e3;            % orbit altitude (m)
inc = (20:2:46)';      % incidence angles (deg)
Re = 6371e3;
% slant range from altitude and incidence (spherical Earth)
look = asin(Re / (Re + H) * sind(inc));
R = (Re + H) * cos(look) - sqrt(Re^2 - ((Re + H) * sin(look)).^2);
Pavg = Pt * p.Tp * p.PRF;
nesz = 2 * (4*pi)^3 .* R.^3 * p.V * k * T0 * F * p.Br * L .* sind(inc) ./ (Pavg * G^2 * p.lambda^3 * p.c);
T = [inc, R / 1e3, 10 * log10(nesz)];
end
