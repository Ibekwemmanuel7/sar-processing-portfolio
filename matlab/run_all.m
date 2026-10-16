%RUN_ALL  Simulate, focus, measure and verify against the requirement and the analytic model.
%   Two processing cases: unweighted (alpha = 1.0) and Hamming-weighted (alpha = 0.75, weighted reference case).
%   Writes results/irf_metrics.csv, results/nesz.csv and results/verification.txt.
p0 = sar_params();
pf = {'FAIL', 'PASS'};
if ~exist('results', 'dir'), mkdir('results'); end
fprintf('C-band stripmap: lambda %.4f m, PRF %.0f Hz, Ba %.0f Hz, Ta %.3f s\n', p0.lambda, p0.PRF, p0.Ba, p0.Ta);
tic; [raw, t, eta, p0] = simulate_raw(p0);
fprintf('raw data %d x %d simulated in %.1f s\n', size(raw, 1), size(raw, 2), toc);
fid = fopen('results/irf_metrics.csv', 'w');
fprintf(fid, 'alpha,range_filter,target,az_m,rg_m,irw_rg_m,irw_az_m,pslr_rg_db,pslr_az_db,islr_rg_db,islr_az_db,loc_err_rg_m,loc_err_az_m\n');
vf = fopen('results/verification.txt', 'w');
allpass = true;
cases = {1.0, 'replica'; 0.75, 'posp'; 0.75, 'replica'};
for c = 1:size(cases, 1)
    alpha = cases{c, 1};
    p = p0; p.win_alpha = alpha; p.range_filter = cases{c, 2};
    tic; [img, info] = rda_focus(raw, t, p);
    mdl = irf_model(p);
    fprintf('\n=== alpha %.2f, %s range filter: focused in %.1f s, max RCM %.2f range cells ===\n', alpha, p.range_filter, toc, info.max_rcm_cells);
    fprintf('model : IRW rg %.2f az %.2f m | PSLR rg %.1f az %.1f dB | ISLR rg %.1f az %.1f dB\n', ...
        mdl.irw_rg_m, mdl.irw_az_m, mdl.pslr_rg_db, mdl.pslr_az_db, mdl.islr_rg_db, mdl.islr_az_db);
    fprintf(vf, 'alpha %.2f model: IRW rg %.3f az %.3f m, PSLR rg %.2f az %.2f dB, ISLR rg %.2f az %.2f dB\n', ...
        alpha, mdl.irw_rg_m, mdl.irw_az_m, mdl.pslr_rg_db, mdl.pslr_az_db, mdl.islr_rg_db, mdl.islr_az_db);
    for k = 1:size(p.targets, 1)
        m = measure_irf(img, p, t, eta, p.targets(k, 1), p.targets(k, 2));
        fprintf(fid, '%.2f,%s,%d,%g,%g,%.3f,%.3f,%.2f,%.2f,%.2f,%.2f,%.3f,%.3f\n', alpha, p.range_filter, k, p.targets(k, 1), p.targets(k, 2), ...
            m.irw_rg_m, m.irw_az_m, m.pslr_rg_db, m.pslr_az_db, m.islr_rg_db, m.islr_az_db, m.loc_err_rg_m, m.loc_err_az_m);
        line = sprintf('T%d    : IRW rg %.2f az %.2f m | PSLR rg %.1f az %.1f dB | ISLR rg %.1f az %.1f dB | loc err %.2f / %.2f m', ...
            k, m.irw_rg_m, m.irw_az_m, m.pslr_rg_db, m.pslr_az_db, m.islr_rg_db, m.islr_az_db, m.loc_err_rg_m, m.loc_err_az_m);
        if alpha == p0.win_alpha && strcmp(p.range_filter, 'replica')
            s = p.spec;
            checks = [m.irw_rg_m <= s.irw_rg_m, m.irw_az_m <= s.irw_az_m, ...
                      max(m.pslr_rg_db, m.pslr_az_db) <= s.pslr_db, max(m.islr_rg_db, m.islr_az_db) <= s.islr_db, ...
                      max(abs([m.loc_err_rg_m m.loc_err_az_m])) <= s.loc_err_m, ...
                      abs(m.irw_rg_m / mdl.irw_rg_m - 1) <= s.model_tol, abs(m.irw_az_m / mdl.irw_az_m - 1) <= s.model_tol];
            ok = all(checks); allpass = allpass && ok;
            line = [line ' | ' pf{ok + 1}];
        end
        fprintf('%s\n', line); fprintf(vf, 'alpha %.2f %s: %s\n', alpha, p.range_filter, line);
    end
    if alpha == p0.win_alpha && strcmp(p.range_filter, 'replica'), save('-v7', 'results/focused.mat', 'img', 't', 'eta'); end
end
fclose(fid);
N = nesz_model(p0);
dlmwrite('results/nesz.csv', N, 'precision', '%.3f');
nesz_ok = max(N(:, 3)) <= p0.spec.nesz_db; allpass = allpass && nesz_ok;
fprintf('\nNESZ across swath (%d to %d deg incidence): %.1f to %.1f dB, requirement <= %.0f dB: %s\n', ...
    N(1, 1), N(end, 1), min(N(:, 3)), max(N(:, 3)), p0.spec.nesz_db, pf{nesz_ok + 1});
fprintf(vf, 'NESZ %.2f to %.2f dB (requirement <= %.1f dB): %s\n', min(N(:, 3)), max(N(:, 3)), p0.spec.nesz_db, pf{nesz_ok + 1});
ov = {'REQUIREMENTS NOT MET', 'ALL REQUIREMENTS MET (weighted case)'};
fprintf('\nOVERALL: %s\n', ov{allpass + 1}); fprintf(vf, 'OVERALL: %s\n', ov{allpass + 1});
fclose(vf);
