function mask = detect_dark(I, k_db, win_bg, win_sm, min_px)
%DETECT_DARK  Adaptive-threshold dark-spot detector for oil slicks in SAR intensity.
%   A pixel is dark if its local mean (win_sm x win_sm) is k_db below the background mean
%   (win_bg x win_bg). Connected regions smaller than min_px are removed.
if nargin < 2, k_db = 3; end
if nargin < 3, win_bg = 101; end
if nargin < 4, win_sm = 7; end
if nargin < 5, min_px = 200; end
bg = boxmean_(I, win_bg); sm = boxmean_(I, win_sm);
mask = sm < bg * 10^(-k_db / 10);
mask = remove_small_(mask, min_px);
end

function m = boxmean_(I, w)
k = ones(w) / w^2;
num = conv2(I, k, 'same'); den = conv2(ones(size(I)), k, 'same');
m = num ./ den;
end

function out = remove_small_(mask, min_px)
% connected components by iterative flood fill (4-connectivity), toolbox-free
[n1, n2] = size(mask); lab = zeros(n1, n2); out = false(n1, n2); cur = 0;
for i = 1:n1*n2
    if mask(i) && lab(i) == 0
        cur = cur + 1; stack = i; lab(i) = cur; members = [];
        while ~isempty(stack)
            q = stack(end); stack(end) = []; members(end+1) = q; %#ok<AGROW>
            [r, c] = ind2sub([n1 n2], q);
            nb = [r-1 c; r+1 c; r c-1; r c+1];
            nb = nb(nb(:,1) >= 1 & nb(:,1) <= n1 & nb(:,2) >= 1 & nb(:,2) <= n2, :);
            for j = 1:size(nb, 1)
                qi = sub2ind([n1 n2], nb(j,1), nb(j,2));
                if mask(qi) && lab(qi) == 0, lab(qi) = cur; stack(end+1) = qi; end %#ok<AGROW>
            end
        end
        if numel(members) >= min_px, out(members) = true; end
    end
end
end
