% Validate symmetry-orbit builders on full FFT R-grid and ISDF sampling grid.
%
% Run from test_Si folder:
%   validate_rgrid_symm_orbits

here = fileparts(mfilename('fullpath'));
cd(here);
cd ../../;
QPstartup;
cd(here);

load('SAVE/config.mat', 'config');
if ~config.ISDF.isisdf
  error('validate_rgrid_symm_orbits:isisdf', ...
    'Set isisdf in test input (./test) before running this validation.');
end

stagePath = 'test_relay_stage.mat';
if ~exist(stagePath, 'file')
  error('validate_rgrid_symm_orbits:stage', ...
    'Missing %s. Run input_driver first.', stagePath);
end

relay.stage_from_db(stagePath);
relay.restore();
isdf.driver([], config);

% Locate current "nn" ISDF slot.
L = isdf.manager('list');
id_nn = [];
for k = 1:numel(L)
  if L(k).assigned && strcmp(char(L(k).desc), 'nn')
    id_nn = L(k).id;
    break
  end
end
if isempty(id_nn)
  error('validate_rgrid_symm_orbits:id', ...
    'No assigned ISDF slot with desc ''nn'' after isdf.driver.');
end

[N1, Nr1, R1, Rrep1, ir2rep1, ir2rot1] = rgrid_symm_orbit();
[N2, Nr2, R2, Rrep2, ir2rep2, ir2rot2] = isdf_r_sampling_symm_orbit(id_nn);
nsym = int32(size(symmetry.get().rot_mtrx_RLU_R, 3));

assert(size(R1, 1) == double(Nr1), 'R1 row size mismatch.');
assert(size(R2, 1) == double(Nr2), 'R2 row size mismatch.');
assert(size(Rrep1, 1) == double(N1), 'Rrep1 row size mismatch.');
assert(size(Rrep2, 1) == double(N2), 'Rrep2 row size mismatch.');
assert(all(ir2rep1 > 0) && all(ir2rep2 > 0), 'ir2rep has uncovered entries.');
assert(all(ir2rot1 > 0) && all(ir2rot2 > 0), 'ir2rot has uncovered entries.');
assert(max(ir2rep1) == N1 && max(ir2rep2) == N2, 'Representative index range mismatch.');
assert(max(ir2rot1) <= nsym && max(ir2rot2) <= nsym, 'Rotation index out of range.');

report = struct();
report.ok = true;
report.nsym = nsym;
report.rgrid = struct('N_r_orbit', N1, 'Nr', Nr1, 'N_rep_rows', int32(size(Rrep1, 1)));
report.isdf_sampling = struct('id', int32(id_nn), 'N_r_orbit', N2, ...
  'Nr', Nr2, 'N_rep_rows', int32(size(Rrep2, 1)));

% Orbit-wise Gram condition number check.
% Gram construction follows isdf_prod convention:
%   G(i,j) = (Psi_i * Psi_j^H) .* (Phi_i * Phi_j^H)
wf_data = wave_functions.get();
isdf_data = isdf.get(id_nn);
[nrange1, nrange2] = isdf_get_nrange(id_nn);

Psi_rgrid = wf_data.c(:, nrange1, 1, 1);
Phi_rgrid = wf_data.c(:, nrange2, 1, 1);
Psi_sampling = isdf_data.coeff_seper(:, nrange1, 1, 1);
Phi_sampling = isdf_data.coeff_seper(:, nrange2, 1, 1);

assert(size(Psi_rgrid, 1) == double(Nr1), 'Psi_rgrid row size mismatch.');
assert(size(Psi_sampling, 1) == double(Nr2), 'Psi_sampling row size mismatch.');

[cond_rgrid, size_rgrid] = local_orbit_gram_condition(Psi_rgrid, Phi_rgrid, ir2rep1, N1);
[cond_sampling, size_sampling] = local_orbit_gram_condition(Psi_sampling, Phi_sampling, ir2rep2, N2);

[max_cond_rgrid, idx_worst_rgrid] = max(cond_rgrid);
[max_cond_sampling, idx_worst_sampling] = max(cond_sampling);
report.rgrid.gram_cond_min = min(cond_rgrid);
report.rgrid.gram_cond_median = median(cond_rgrid);
report.rgrid.gram_cond_max = max_cond_rgrid;
report.rgrid.worst_orbit = int32(idx_worst_rgrid);
report.rgrid.worst_orbit_size = int32(size_rgrid(idx_worst_rgrid));

report.isdf_sampling.gram_cond_min = min(cond_sampling);
report.isdf_sampling.gram_cond_median = median(cond_sampling);
report.isdf_sampling.gram_cond_max = max_cond_sampling;
report.isdf_sampling.worst_orbit = int32(idx_worst_sampling);
report.isdf_sampling.worst_orbit_size = int32(size_sampling(idx_worst_sampling));

log_path = fullfile(here, 'SAVE', 'orbit_condition_log.txt');
local_write_orbit_cond_log(log_path, report, size_rgrid, cond_rgrid, size_sampling, cond_sampling);

fprintf('\n=== validate_rgrid_symm_orbits report ===\n');
fprintf('status: OK\n');
fprintf('nsym: %d\n', double(report.nsym));
fprintf('rgrid_symm_orbit: N_r_orbit=%d, Nr=%d, rep_rows=%d\n', ...
  double(report.rgrid.N_r_orbit), double(report.rgrid.Nr), double(report.rgrid.N_rep_rows));
fprintf('isdf_r_sampling_symm_orbit (id=%d): N_r_orbit=%d, Nr=%d, rep_rows=%d\n', ...
  double(report.isdf_sampling.id), ...
  double(report.isdf_sampling.N_r_orbit), ...
  double(report.isdf_sampling.Nr), ...
  double(report.isdf_sampling.N_rep_rows));
fprintf('rgrid Gram cond per orbit: min=%.3e, median=%.3e, max=%.3e (worst orbit=%d, size=%d)\n', ...
  report.rgrid.gram_cond_min, ...
  report.rgrid.gram_cond_median, ...
  report.rgrid.gram_cond_max, ...
  double(report.rgrid.worst_orbit), ...
  double(report.rgrid.worst_orbit_size));
fprintf('isdf sampling Gram cond per orbit: min=%.3e, median=%.3e, max=%.3e (worst orbit=%d, size=%d)\n', ...
  report.isdf_sampling.gram_cond_min, ...
  report.isdf_sampling.gram_cond_median, ...
  report.isdf_sampling.gram_cond_max, ...
  double(report.isdf_sampling.worst_orbit), ...
  double(report.isdf_sampling.worst_orbit_size));
fprintf('orbit condition log saved: %s\n', log_path);
fprintf('=== end report ===\n\n');



function [cond_per_orbit, orbit_sizes] = local_orbit_gram_condition(Psi_all, Phi_all, ir2rep, N_orbit)
  n_orbit = double(N_orbit);
  cond_per_orbit = zeros(n_orbit, 1);
  orbit_sizes = zeros(n_orbit, 1);

  for irep = 1:n_orbit
    idx = find(double(ir2rep) == irep);
    orbit_sizes(irep) = numel(idx);
    if isempty(idx)
      cond_per_orbit(irep) = Inf;
      continue;
    end

    Psi_orbit = Psi_all(idx, :);
    Phi_orbit = Phi_all(idx, :);
    norb = size(Psi_orbit, 1);
    G = zeros(norb, norb);
    G = isdf_prod(Psi_orbit, Psi_orbit, Phi_orbit, Phi_orbit);
    G = (G + G') / 2;

    s = svd(double(G));
    if isempty(s) || s(1) == 0
      cond_per_orbit(irep) = Inf;
    else
      cond_per_orbit(irep) = s(1) / max(s(end), eps(s(1)));
    end
    if cond_per_orbit(irep) > 1e+6
      ;
    end
  end
end

function local_write_orbit_cond_log(log_path, report, size_rgrid, cond_rgrid, size_sampling, cond_sampling)
  [log_dir, ~, ~] = fileparts(log_path);
  if ~exist(log_dir, 'dir')
    mkdir(log_dir);
  end

  fid = fopen(log_path, 'w');
  if fid < 0
    error('validate_rgrid_symm_orbits:LogOpenFailed', ...
      'Cannot open log file for writing: %s', log_path);
  end
  cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

  fprintf(fid, 'validate_rgrid_symm_orbits orbit-condition log\n');
  fprintf(fid, 'status: OK\n');
  fprintf(fid, 'nsym: %d\n', double(report.nsym));
  fprintf(fid, '\n');

  fprintf(fid, '[rgrid_symm_orbit]\n');
  fprintf(fid, 'N_r_orbit=%d, Nr=%d\n', ...
    double(report.rgrid.N_r_orbit), double(report.rgrid.Nr));
  fprintf(fid, '%-14s %-12s %-18s\n', 'orbit_number', 'orbit_size', 'condition_number');
  for i = 1:numel(cond_rgrid)
    fprintf(fid, '%-14d %-12d %-18.10e\n', i, size_rgrid(i), cond_rgrid(i));
  end
  fprintf(fid, '\n');

  fprintf(fid, '[isdf_r_sampling_symm_orbit]\n');
  fprintf(fid, 'id=%d, N_r_orbit=%d, Nr=%d\n', ...
    double(report.isdf_sampling.id), ...
    double(report.isdf_sampling.N_r_orbit), ...
    double(report.isdf_sampling.Nr));
  fprintf(fid, '%-14s %-12s %-18s\n', 'orbit_number', 'orbit_size', 'condition_number');
  for i = 1:numel(cond_sampling)
    fprintf(fid, '%-14d %-12d %-18.10e\n', i, size_sampling(i), cond_sampling(i));
  end
end

