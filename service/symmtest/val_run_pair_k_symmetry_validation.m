% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/19 ZZ

function report = val_run_pair_k_symmetry_validation(ctx, opts)
% Composite k = (k1, k2) on BZ mesh: pair-product equivariance under S.
% Compare rotate(conj(psi1).*psi2) with conj(psi1').*psi2' for ki' = S*ki.

if nargin < 2
  opts = struct();
end
if ~isfield(opts, 'tol')
  opts.tol = 1e-5;
end
if ~isfield(opts, 'max_bands')
  opts.max_bands = min(4, double(ctx.wf_data.nb));
end
if ~isfield(opts, 'max_pairs')
  opts.max_pairs = 24;
end
if ~isfield(opts, 'pair_k_publish_base')
  opts.pair_k_publish_base = true;
end

tol = opts.tol;
k_data = ctx.k_data;
symm_data = ctx.symm_data;
fft_data = ctx.fft_data;
wf_data = ctx.wf_data;
system_data = system.manager('get');

kptbz = double(k_data.kptbz_RLU);
kptbz_Cart = double(k_data.kptbz_Cart);
nkbz = size(kptbz, 1);
nsym = size(symm_data.rot_mtrx_RLU_G, 3);
nspin = double(wf_data.nspin);
nb = min(double(wf_data.nb), opts.max_bands);

pairs = zeros(0, 2, 'int32');
for ik1 = 1:nkbz
  for ik2 = ik1:nkbz
    pairs(end + 1, :) = int32([ik1, ik2]); %#ok<AGROW>
    if size(pairs, 1) >= opts.max_pairs
      break;
    end
  end
  if size(pairs, 1) >= opts.max_pairs
    break;
  end
end

d_err = [];
fails = {};
debug_records = {};

for ispin = 1:nspin
  for ib = 1:nb
    for ip = 1:size(pairs, 1)
      ik1 = pairs(ip, 1);
      ik2 = pairs(ip, 2);

      wf1 = get_wf_at_bz(ib, ik1, ispin, k_data);
      wf2 = get_wf_at_bz(ib, ik2, ispin, k_data);

      % Bamp-style real-space product used in SCATTER_Bamp.
      % phi(r) = conj(u_{k1,b}(r)) .* u_{k2,b}(r)
      prod_field = conj(wf1) .* wf2;

      for isym = 1:nsym
        S_rlu = double(symm_data.rot_mtrx_RLU_G(:, :, isym));
        S_cart = double(symm_data.rot_mtrx_Cart(:, :, isym));
        % Sk in both coordinate systems (for diagnostics and G0 extraction).
        % In RLU: Sk = k * S, then k' = Sk + G0, where G0 is integer (RLU).
        Sk1_RLU = double(kptbz(ik1, :)) * S_rlu;
        Sk2_RLU = double(kptbz(ik2, :)) * S_rlu;
        Sk1_Cart = double(kptbz_Cart(ik1, :)) * S_cart;
        Sk2_Cart = double(kptbz_Cart(ik2, :)) * S_cart;

        ik1p = map_kbz_by_symm(ik1, isym, kptbz, symm_data.rot_mtrx_RLU_G, tol);
        ik2p = map_kbz_by_symm(ik2, isym, kptbz, symm_data.rot_mtrx_RLU_G, tol);
        % G0_i = round(k_i' - S*k_i) in RLU, used only for reporting here.
        if ik1p > 0
          G0_1 = round(double(kptbz(ik1p, :)) - Sk1_RLU);
        else
          G0_1 = [];
        end
        if ik2p > 0
          G0_2 = round(double(kptbz(ik2p, :)) - Sk2_RLU);
        else
          G0_2 = [];
        end

        if ik1p < 1 || ik2p < 1
          msg = sprintf('pair map missing: pair=(%d,%d), isym=%d', ik1, ik2, isym);
          fails{end + 1} = msg; %#ok<AGROW>
          rec = struct( ...
            'kind', 'map_missing', ...
            'fail_message', msg, ...
            'ik1', int32(ik1), ...
            'ik2', int32(ik2), ...
            'ik1p', int32(ik1p), ...
            'ik2p', int32(ik2p), ...
            'ib', int32(ib), ...
            'ispin', int32(ispin), ...
            'isym', int32(isym), ...
            'rel_err', nan, ...
            'k1_Cart', kptbz_Cart(ik1, :), ...
            'k2_Cart', kptbz_Cart(ik2, :), ...
            'k1p_Cart', [], ...
            'k2p_Cart', [], ...
            'Sk1_Cart', Sk1_Cart, ...
            'Sk2_Cart', Sk2_Cart, ...
            'G0_1_RLU', G0_1, ...
            'G0_2_RLU', G0_2, ...
            'S_matrix', S_cart, ...
            'ib_list_1', int32([]), ...
            'ib_list_2', int32([]), ...
            'lhs', [], ...
            'Vpair', []);
          debug_records{end + 1} = rec; %#ok<AGROW>
          if opts.pair_k_publish_base
            pair_k_print_and_publish_debug(rec);
          end
          continue;
        end

        % lhs = S[ conj(u_{k1,b}) .* u_{k2,b} ]
        lhs = apply_symm_on_field(prod_field, isym, symm_data, fft_data);

        ik1ibz = double(k_data.bz2ibz(ik1));
        ik2ibz = double(k_data.bz2ibz(ik2));
        ib_list_1 = get_degen_band_indices(ib, ik1ibz, system_data, wf_data);
        ib_list_2 = get_degen_band_indices(ib, ik2ibz, system_data, wf_data);
        % Vpair spans products from degen spaces at mapped points:
        % span{ conj(u_{k1p,b1}) .* u_{k2p,b2} | b1 in D1, b2 in D2 }.
        Vpair = build_pair_subspace(ib_list_1, ib_list_2, ik1p, ik2p, ispin, k_data);

        % Relative distance from lhs to pair subspace:
        % rel_err = ||lhs - P_V(lhs)|| / ||lhs||, P_V = Vpair*Vpair'.
        rel_err = subspace_projection_rel_error(lhs, Vpair);
        d_err(end + 1, 1) = rel_err; %#ok<AGROW>

        if rel_err > 50 * tol
          msg = sprintf('pair-product fail: ib=%d pair=(%d,%d) isym=%d err=%.3e', ib, ik1, ik2, isym, rel_err);
          fails{end + 1} = msg; %#ok<AGROW>
          rec = struct( ...
            'kind', 'pair_subspace_mismatch', ...
            'fail_message', msg, ...
            'ik1', int32(ik1), ...
            'ik2', int32(ik2), ...
            'ik1p', int32(ik1p), ...
            'ik2p', int32(ik2p), ...
            'ib', int32(ib), ...
            'ispin', int32(ispin), ...
            'isym', int32(isym), ...
            'rel_err', rel_err, ...
            'k1_Cart', kptbz_Cart(ik1, :), ...
            'k2_Cart', kptbz_Cart(ik2, :), ...
            'k1p_Cart', kptbz_Cart(ik1p, :), ...
            'k2p_Cart', kptbz_Cart(ik2p, :), ...
            'Sk1_Cart', Sk1_Cart, ...
            'Sk2_Cart', Sk2_Cart, ...
            'G0_1_RLU', G0_1, ...
            'G0_2_RLU', G0_2, ...
            'S_matrix', S_cart, ...
            'ib_list_1', ib_list_1, ...
            'ib_list_2', ib_list_2, ...
            'lhs', lhs, ...
            'Vpair', Vpair);
          debug_records{end + 1} = rec; %#ok<AGROW>
          if opts.pair_k_publish_base
            pair_k_print_and_publish_debug(rec);
          end
        end
      end
    end
  end
end

report = build_report_struct('pair_product', d_err, fails, tol);
report.n_pairs_tested = int32(size(pairs, 1));

if opts.pair_k_publish_base && ~isempty(debug_records)
  assignin('base', 'pair_k_debug_all', debug_records);
end
end

function pair_k_print_and_publish_debug(rec)
% Print pair-k failure details and publish latest record to base workspace.

fprintf('\n=== pair-k validation failure ===\n');
fprintf('kind:          %s\n', rec.kind);
fprintf('fail_message:  %s\n', rec.fail_message);
fprintf('ik1=%d ik2=%d -> ik1p=%d ik2p=%d  ib=%d ispin=%d isym=%d  rel_err=%.6g\n', ...
  rec.ik1, rec.ik2, rec.ik1p, rec.ik2p, rec.ib, rec.ispin, rec.isym, rec.rel_err);
fprintf('k1 (Cart):  %s\n', mat2str(rec.k1_Cart, 8));
fprintf('k2 (Cart):  %s\n', mat2str(rec.k2_Cart, 8));
if ~isempty(rec.k1p_Cart)
  fprintf('k1p (Cart): %s\n', mat2str(rec.k1p_Cart, 8));
else
  fprintf('k1p (Cart): (n/a, map missing)\n');
end
if ~isempty(rec.k2p_Cart)
  fprintf('k2p (Cart): %s\n', mat2str(rec.k2p_Cart, 8));
else
  fprintf('k2p (Cart): (n/a, map missing)\n');
end
fprintf('Sk1 (Cart): %s\n', mat2str(rec.Sk1_Cart, 8));
fprintf('Sk2 (Cart): %s\n', mat2str(rec.Sk2_Cart, 8));
if ~isempty(rec.G0_1_RLU)
  fprintf('G0_1 (RLU): [%s]\n', num2str(double(rec.G0_1_RLU(:)')));
else
  fprintf('G0_1 (RLU): (n/a, map missing)\n');
end
if ~isempty(rec.G0_2_RLU)
  fprintf('G0_2 (RLU): [%s]\n', num2str(double(rec.G0_2_RLU(:)')));
else
  fprintf('G0_2 (RLU): (n/a, map missing)\n');
end
fprintf('S_matrix (3x3, rot_mtrx_Cart(:,:,isym)):\n');
disp(rec.S_matrix);
if ~isempty(rec.lhs)
  fprintf('lhs:   size %s  ||.||_2 = %.6g\n', mat2str(size(rec.lhs)), norm(rec.lhs(:)));
  fprintf('Vpair: size %s\n', mat2str(size(rec.Vpair)));
  fprintf('ib_list_1: [%s]\n', num2str(double(rec.ib_list_1(:)')));
  fprintf('ib_list_2: [%s]\n', num2str(double(rec.ib_list_2(:)')));
  proj = rec.Vpair * (rec.Vpair' * rec.lhs(:));
  fprintf('||lhs - Vpair*Vpair''*lhs||_2: %.6g\n', norm(rec.lhs(:) - proj));
else
  fprintf('lhs / Vpair: (empty for this failure kind)\n');
end
fprintf('==================================\n\n');

assignin('base', 'pair_k_debug_latest', rec);
end

function ib_list = get_degen_band_indices(ib, ikibz, system_data, wf_data)
% Return band indices that belong to ib's degeneracy segment at ikibz.
% If ib is in [ib_start, ib_start+n_degen-1], return that full segment.
nb_total = double(wf_data.nb);
first_arr = system_data.first_index_in_degeneracy{ikibz};
num_arr = system_data.num_index_in_degeneracy{ikibz};

ib_start = ib;
n_degen = 1;
for iseg = 1:numel(first_arr)
  seg_start = double(first_arr(iseg));
  seg_len = double(num_arr(iseg));
  if ib >= seg_start && ib < seg_start + seg_len
    ib_start = seg_start;
    n_degen = seg_len;
    break;
  end
end

ib_end = min(ib_start + n_degen - 1, nb_total);
ib_list = int32(ib_start:ib_end);
end

function Vpair = build_pair_subspace(ib_list_1, ib_list_2, ik1p, ik2p, ispin, k_data)
% Build pair-product candidate matrix and orthonormalize its column space.
% Raw columns: v_{b1,b2}(r) = conj(u_{k1p,b1}(r)) .* u_{k2p,b2}(r).
cols = cell(0, 1);
for i1 = 1:numel(ib_list_1)
  wf1p = get_wf_at_bz(double(ib_list_1(i1)), ik1p, ispin, k_data);
  for i2 = 1:numel(ib_list_2)
    wf2p = get_wf_at_bz(double(ib_list_2(i2)), ik2p, ispin, k_data);
    cols{end + 1, 1} = conj(wf1p(:)) .* wf2p(:); %#ok<AGROW>
  end
end

if isempty(cols)
  Vpair = zeros(0, 0);
  return;
end

% Explicitly concatenate candidate basis vectors as columns.
Vraw = double(cat(2, cols{:}));
% Vpair has orthonormal columns; Vpair*Vpair' is the orthogonal projector.
Vpair = orth(Vraw);
end

function rel_err = subspace_projection_rel_error(lhs, V)
% Relative residual of lhs outside the column space of V.
% rel_err = ||lhs - V*(V'*lhs)|| / ||lhs||.
lhs = lhs(:);
nrm = norm(lhs);
if nrm < eps
  rel_err = 0;
  return;
end
if isempty(V)
  rel_err = 1;
  return;
end
rel_err = norm(lhs - V * (V' * lhs)) / nrm;
end
