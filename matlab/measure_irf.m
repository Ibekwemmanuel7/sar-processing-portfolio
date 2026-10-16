function m = measure_irf(img, p, t, eta, xa, dr)
%MEASURE_IRF  Impulse-response metrics around the target at azimuth xa (m), range offset dr (m).
%   Upsamples a 32 x 32 patch by 16 (FFT zero padding), then measures -3 dB width, PSLR and ISLR
%   along range and azimuth cuts through the peak, and the peak location error.
c = p.c; dx_r = c / (2 * p.fs); dx_a = p.V / p.PRF;
R = t * c / 2;
[~, ir] = min(abs(R - (p.R0 + dr)));
[~, ia] = min(abs(p.V * eta - xa));
[~, k] = max(abs(img(ia-20:ia+20, ir-20:ir+20)(:)));        % local peak
[ka, kr] = ind2sub([41 41], k); ia = ia - 21 + ka; ir = ir - 21 + kr;
h = 16; n = 32;
patch = img(ia-n/2:ia+n/2-1, ir-n/2:ir+n/2-1);
P = fftshift(fft2(patch));
Pz = zeros(n*h, n*h); o = (n*h - n)/2;
Pz(o+1:o+n, o+1:o+n) = P;
up = abs(ifft2(ifftshift(Pz))).^2;
[~, k] = max(up(:)); [pa, pr] = ind2sub(size(up), k);
cut_r = up(pa, :) / up(pa, pr); cut_a = up(:, pr)' / up(pa, pr);
[m.irw_rg_m, m.pslr_rg_db, m.islr_rg_db] = cut_metrics(cut_r, pr, dx_r / h);
[m.irw_az_m, m.pslr_az_db, m.islr_az_db] = cut_metrics(cut_a, pa, dx_a / h);
% location error: peak position vs true position
pk_r = R(ir - n/2) + (pr - 1) * dx_r / h;
pk_a = p.V * eta(ia - n/2) + (pa - 1) * dx_a / h;
m.loc_err_rg_m = pk_r - (p.R0 + dr);
m.loc_err_az_m = pk_a - xa;
m.peak = sqrt(up(pa, pr));
end

