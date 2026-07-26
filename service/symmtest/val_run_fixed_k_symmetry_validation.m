function report = val_run_fixed_k_symmetry_validation(ctx, opts)
% Single-BZ-index k: stabilizer check (exists S with S*k = k mod G*).
% Compare apply_symm(psi(k)) with psi at the mesh index for S*k.

if nargin < 2
  opts = struct();
end
if ~isfield(opts, 'tol')
  opts.tol = 1e-5;
end
if ~isfield(opts, 'max_bands')
  opts.max_bands = min(4, double(ctx.wf_data.nb));
end
if ~isfield(opts, 'fixed_k_skip_identity')
  opts.fixed_k_skip_identity = true;
end
if ~isfield(opts, 'save_fixed_k_failures')
  opts.save_fixed_k_failures = true;
end
if ~isfield(opts, 'fixed_k_debug_dir')
  opts.fixed_k_debug_dir = '';
end
if ~isfield(opts, 'fixed_k_publish_base')
  opts.fixed_k_publish_base = true;
end

tol = opts.tol;
k_data = ctx.k_data;
symm_data = ctx.symm_data;
fft_data = ctx.fft_data;
wf_data = ctx.wf_data;
system_data = system.manager('get');
rotG = symm_data.rot_mtrx_RLU_G;

kptbz = double(k_data.kptbz_RLU);
kptbz_Cart = double(k_data.kptbz_Cart);
nkbz = size(kptbz, 1);
nsym = size(rotG, 3);
nspin = double(wf_data.nspin);
nb = min(double(wf_data.nb), opts.max_bands);

d_err = [];
fails = {};
debug_records = {};
n_stab_checks = int32(0);

for ik = 1:nkbz
  for isym = 1:nsym
    if opts.fixed_k_skip_identity && isym == 1
      continue;
    end
    if ~k_point_fixed_by_symmetry(ik, isym, kptbz, rotG, tol)
      continue;
    end

    ik_map = map_kbz_by_symm(ik, isym, kptbz, rotG, tol);
    if ik_map < 1
      S_mat = double(symm_data.rot_mtrx_Cart(:, :, isym));
      msg = sprintf('fixed-k map missing: ik=%d isym=%d', ik, isym);
      fails{end + 1} = msg; %#ok<AGROW>
      rec = struct( ...
        'kind', 'map_missing', ...
        'fail_message', msg, ...
        'ik', int32(ik), ...
        'ik_map', int32(-1), ...
        'ib', int32(-1), ...
        'ispin', int32(-1), ...
        'isym', int32(isym), ...
        'rel_err', nan, ...
        'k_Cart', kptbz_Cart(ik, :), ...
        'k_map_Cart', [], ...
        'S_matrix', S_mat, ...
        'wf_k', [], ...
        'lhs', [], ...
        'V', []);
      debug_records{end + 1} = rec; %#ok<AGROW>
      if opts.fixed_k_publish_base
        fixed_k_print_and_publish_debug(rec);
      end
      continue;
    end

    ikibz = double(k_data.bz2ibz(ik));
    [k_shift_G0, has_k_shift] = get_kmap_shift_g0(ik, ik_map, isym, kptbz, rotG, tol);
    for ispin = 1:nspin
      for ib = 1:nb
        wf_k = get_wf_at_bz(ib, ik, ispin, k_data);
        lhs = apply_symm_on_field(wf_k, isym, symm_data, fft_data);
        if has_k_shift
          lhs = apply_kshift_phase_back(lhs, k_shift_G0, fft_data);
        end
        [V, ib_degen_range] = build_degen_subspace(ib, ikibz, ispin, ik_map, system_data, k_data, wf_data);

        

        rel_err = subspace_projection_rel_error(lhs, V);
        d_err(end + 1, 1) = rel_err; %#ok<AGROW>
        n_stab_checks = n_stab_checks + 1;

        if rel_err > 20 * tol
          S_mat = double(symm_data.rot_mtrx_Cart(:, :, isym));
          msg = sprintf( ...
            'fixed-k fail: ik=%d ik_map=%d ib=%d ispin=%d isym=%d err=%.3e', ...
            ik, ik_map, ib, ispin, isym, rel_err);
          fails{end + 1} = msg; %#ok<AGROW>
          rec = struct( ...
            'kind', 'wf_mismatch', ...
            'fail_message', msg, ...
            'ik', int32(ik), ...
            'ik_map', int32(ik_map), ...
            'ib', int32(ib), ...
            'ispin', int32(ispin), ...
            'isym', int32(isym), ...
            'k_shift_G0', k_shift_G0, ...
            'degen_band_indices', ib_degen_range, ...
            'rel_err', rel_err, ...
            'k_Cart', kptbz_Cart(ik, :), ...
            'k_map_Cart', kptbz_Cart(ik_map, :), ...
            'S_matrix', S_mat, ...
            'wf_k', wf_k, ...
            'lhs', lhs, ...
            'V', V);
          debug_records{end + 1} = rec; %#ok<AGROW>
          if opts.fixed_k_publish_base
            fixed_k_print_and_publish_debug(rec);
          end
        end
      end
    end
  end
end

report = build_report_struct('fixed_k', d_err, fails, tol);
report.nkbz = int32(nkbz);
report.n_stabilizer_checks = n_stab_checks;
report.debug_matfile = '';
report.debug_dir_used = '';

if opts.save_fixed_k_failures && ~isempty(debug_records)
  dbgdir = opts.fixed_k_debug_dir;
  if isempty(strtrim(dbgdir))
    dbgdir = fullfile(fileparts(mfilename('fullpath')), 'fixed_k_failure_dumps');
  end
  if ~exist(dbgdir, 'dir')
    mkdir(dbgdir);
  end
  stamp = datestr(now, 'yyyymmdd_HHMMSS');
  matpath = fullfile(dbgdir, sprintf('fixed_k_failures_%s.mat', stamp));
  fixed_k_debug_records = debug_records; %#ok<NASGU>
  fixed_k_fails = fails; %#ok<NASGU>
  fixed_k_tol = tol; %#ok<NASGU>
  fixed_k_nkbz = int32(nkbz); %#ok<NASGU>
  save(matpath, 'fixed_k_debug_records', 'fixed_k_fails', 'fixed_k_tol', 'fixed_k_nkbz', '-v7.3');
  report.debug_matfile = matpath;
  report.debug_dir_used = dbgdir;
end

if opts.fixed_k_publish_base && ~isempty(debug_records)
  assignin('base', 'fixed_k_debug_all', debug_records);
end
end

function fixed_k_print_and_publish_debug(rec)
% Print k, S, and scalar metadata; publish full struct to base workspace (MATLAB Workspace / Command Window).

fprintf('\n=== fixed-k validation failure ===\n');
fprintf('kind:          %s\n', rec.kind);
fprintf('fail_message:  %s\n', rec.fail_message);
fprintf('ik=%d  ik_map=%d  ib=%d  ispin=%d  isym=%d  rel_err=%.6g\n', ...
  rec.ik, rec.ik_map, rec.ib, rec.ispin, rec.isym, rec.rel_err);
fprintf('k (Cart):     %s\n', mat2str(rec.k_Cart, 8));
if ~isempty(rec.k_map_Cart)
  fprintf('k_map (Cart): %s\n', mat2str(rec.k_map_Cart, 8));
else
  fprintf('k_map (Cart): (n/a, map missing)\n');
end
fprintf('S_matrix (3x3, rot_mtrx_Cart(:,:,isym)):\n');
disp(rec.S_matrix);
if ~isempty(rec.wf_k)
  fprintf('wf_k:  size %s  ||.||_2 = %.6g\n', mat2str(size(rec.wf_k)), norm(rec.wf_k(:)));
  fprintf('lhs:   size %s  ||.||_2 = %.6g\n', mat2str(size(rec.lhs)), norm(rec.lhs(:)));
  if isfield(rec, 'k_shift_G0') && ~isempty(rec.k_shift_G0)
    fprintf('G0 shift (RLU): [%s]\n', num2str(double(rec.k_shift_G0(:)')));
  end
  fprintf('V:     size %s  (degen subspace, %d column(s))\n', mat2str(size(rec.V)), size(rec.V, 2));
  fprintf('degen bands: [%s]\n', num2str(double(rec.degen_band_indices(:)')));
  proj = rec.V * (rec.V' * rec.lhs(:));
  fprintf('||lhs - V*V''*lhs||_2: %.6g\n', norm(rec.lhs(:) - proj));
else
  fprintf('wf_k / lhs / V: (empty for this failure kind)\n');
end
fprintf('==================================\n\n');

assignin('base', 'fixed_k_debug_latest', rec);
end

function [k_shift_G0, has_k_shift] = get_kmap_shift_g0(ik, ik_map, isym, kptbz, rotG, tol)
% Compute reciprocal shift G0 such that k_map = S*k + G0 (fractional RLUs).
S = double(rotG(:, :, isym));
Sk = double(kptbz(ik, :)) * S;
diff = double(kptbz(ik_map, :)) - Sk;
k_shift_G0 = round(diff);
has_k_shift = norm(k_shift_G0) > tol;
end

function out = apply_kshift_phase_back(lhs, k_shift_G0, fft_data)
% Remove exp(i G0 r) phase from lhs on FFT real-space grid.
grid = double(fft_data.fftgrid(:)');
r_rlu = double(fft_data.Rgrid_RLU);
phase_arg = 2 * pi * (r_rlu(:, 1) * (k_shift_G0(1) / grid(1)) + ...
  r_rlu(:, 2) * (k_shift_G0(2) / grid(2)) + ...
  r_rlu(:, 3) * (k_shift_G0(3) / grid(3)));
phase_back = exp(-1i * phase_arg);
out = lhs(:) .* phase_back;
end

function [V, ib_degen_range] = build_degen_subspace(ib, ikibz, ispin, ik_map, system_data, k_data, wf_data)
% Collect wave functions at ik_map belonging to the same degeneracy group
% as band ib at IBZ k-index ikibz, then return orthonormal basis V and the
% corresponding band index range ib_degen_range.
nb_total = double(wf_data.nb);

first_arr = system_data.first_index_in_degeneracy{ikibz};
num_arr   = system_data.num_index_in_degeneracy{ikibz};

ib_start = ib;
n_degen  = 1;
for iseg = 1:numel(first_arr)
  seg_start = double(first_arr(iseg));
  seg_len   = double(num_arr(iseg));
  if ib >= seg_start && ib < seg_start + seg_len
    ib_start = seg_start;
    n_degen  = seg_len;
    break;
  end
end

ib_end = min(ib_start + n_degen - 1, nb_total);
ib_degen_range = int32(ib_start : ib_end);

wf_cols = {};
for ib_d = ib_start : ib_end
  wf_cols{end + 1} = get_wf_at_bz(ib_d, ik_map, ispin, k_data); %#ok<AGROW>
end

V = orth(double(cell2mat(wf_cols)));
end

function rel_err = subspace_projection_rel_error(lhs, V)
% Relative residual of lhs outside the column space of V.
% rel_err = ||lhs - V*(V'*lhs)|| / ||lhs||
lhs = lhs(:);
nrm = norm(lhs);
if nrm < eps
  rel_err = 0;
  return;
end
rel_err = norm(lhs - V * (V' * lhs)) / nrm;
end
