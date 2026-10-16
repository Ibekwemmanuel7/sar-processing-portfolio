function [irw, pslr, islr] = cut_metrics(cut, k0, dx)
%CUT_METRICS  -3 dB width, peak sidelobe ratio and integrated sidelobe ratio of a 1-D power cut.
%   cut: normalised power (peak = 1) at index k0; dx: sample spacing (m).
db = 10 * log10(cut + eps);
% -3 dB width with linear interpolation
l = k0; while l > 1 && db(l) > -3, l = l - 1; end
r = k0; while r < numel(db) && db(r) > -3, r = r + 1; end
xl = l + (-3 - db(l)) / (db(l+1) - db(l));
xr = r - 1 + (-3 - db(r-1)) / (db(r) - db(r-1));
irw = (xr - xl) * dx;
% first nulls around the main lobe
a = k0; while a > 2 && cut(a-1) < cut(a), a = a - 1; end
b = k0; while b < numel(cut) - 1 && cut(b+1) < cut(b), b = b + 1; end
side = [cut(1:a-1) cut(b+1:end)];
pslr = 10 * log10(max(side));
islr = 10 * log10(sum(side) / sum(cut(a:b)));
end
