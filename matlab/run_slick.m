%RUN_SLICK  Verify the dark-spot detector against simulated truth across slick contrast and looks.
%   Reports probability of detection (slick pixels found), false-alarm rate on open sea,
%   and whether the low-wind look-alike is flagged. Writes results/slick_sweep.csv.
if ~exist('results', 'dir'), mkdir('results'); end
fid = fopen('results/slick_sweep.csv', 'w');
fprintf(fid, 'looks,contrast_db,pod,far_open_sea,lookalike_flagged_frac\n');
fprintf('looks  contrast  POD    FAR(open sea)  look-alike flagged\n');
for L = [1 4]
    for c = [2 3 4 6 8 10]
        s.looks = L; s.contrast_db = c; s.seed = 3;
        [I, truth, la, s] = slick_scene(s);
        m = detect_dark(I, 3, 101, 7, 200);
        sea = ~truth & ~la;
        pod = sum(m(:) & truth(:)) / sum(truth(:));
        far = sum(m(:) & sea(:)) / sum(sea(:));
        laf = sum(m(:) & la(:)) / sum(la(:));
        fprintf(fid, '%d,%g,%.4f,%.5f,%.4f\n', L, c, pod, far, laf);
        fprintf('%4d  %6g dB  %.3f  %.5f        %.3f\n', L, c, pod, far, laf);
    end
end
fclose(fid);
