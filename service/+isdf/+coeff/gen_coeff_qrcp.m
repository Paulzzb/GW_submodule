% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/03
%
% QRCP index route ported from kssolvGW .../ISDF/isdf_indices.m (case 'qrcp').

function gen_coeff_qrcp(cfg, id)
%GEN_COEFF_QRCP  Select ISDF sampling points by randomized QRCP and fill slot.
%
%   isdf.coeff.gen_coeff_qrcp(cfg, id)
%
% Algorithm (legacy isdf_indices / exxmethod='qrcp'):
%   1. Build real-space Psi/Phi on the fine FFT grid from nrange1/nrange2 (all BZ k).
%   2. Rank rk = isdf_data.nisdf (set by driver from isdf_ratio * sqrt(N1*N2)).
%   3. Gaussian-sketch Phi/Psi, form product states, QR with column pivoting;
%      take the first rk pivoted row indices as fine-grid sampling points.
%   4. Write R_sampling_RLU / coeff_seper / bundle_struct (interp_scheme='qrcp').
%
% cfg may supply:
%   .seed  - RNG seed (default 0, same as config.ISDF.seed default)

  if nargin < 1 || isempty(cfg)
    cfg = struct();
  end
  if nargin < 2 || isempty(id)
    error('isdf:gen_coeff_qrcp:Id', 'ISDF id is required.');
  end

  fft_data = FFT.get();
  wf_data = wave_functions.get();
  symm_data = symmetry.get();
  k_data = lattice.manager('k', 'get');
  isdf_data = isdf.get(id);

  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    error('isdf:gen_coeff_qrcp:nrange', ...
      'Missing cached nrange in ISDF id=%d. Build it first via isdf.set_nrange(id, config.SYSTEM).', ...
      int32(id));
  end

  nrange1 = double(isdf_data.nrange1(:).');
  nrange2 = double(isdf_data.nrange2(:).');
  Nw = double(wf_data.nc);
  nbz = double(k_data.nbz);

  % --- Build Phi/Psi on the fine grid (same layout as gen_coeff_pseudo) ---
  Psi = zeros(Nw, numel(nrange1) * nbz, 'double');
  Phi = zeros(Nw, numel(nrange2) * nbz, 'double');
  count1 = 1;
  count2 = 1;
  ispin = 1;
  for ikbz = 1:nbz
    ikibz = k_data.bz2ibz(ikbz, 1);
    ikrot = k_data.bz2rot(ikbz, 1);
    for ib1 = 1:numel(nrange1)
      Psi(:, count1) = double(wave_functions.WF_apply_symm([nrange1(ib1), ikibz, ikrot, ispin]));
      count1 = count1 + 1;
    end
    for ib2 = 1:numel(nrange2)
      Phi(:, count2) = double(wave_functions.WF_apply_symm([nrange2(ib2), ikibz, ikrot, ispin]));
      count2 = count2 + 1;
    end
  end

  [m1, n1] = size(Psi);
  [m2, n2] = size(Phi);
  if m1 ~= m2
    error('isdf:gen_coeff_qrcp:Size', 'Psi/Phi row dimensions do not match.');
  end

  rk = double(isdf_data.nisdf);
  if ~(isfinite(rk) && rk >= 1)
    error('isdf:gen_coeff_qrcp:Rank', ...
      'isdf_data.nisdf must be a positive rank before QRCP (id=%d).', int32(id));
  end
  rk = min([rk, m1, n1 * n2]);
  rk = max(1, round(rk));

  seed = 0;
  if isstruct(cfg) && isfield(cfg, 'seed') && ~isempty(cfg.seed)
    seed = double(cfg.seed);
  end
  rng(seed, 'twister');

  % --- Randomized QRCP (legacy isdf_indices case 'qrcp') ---
  rsamp = 1.2 * rk;
  r1 = min(ceil(sqrt((n1 / n2) * rsamp)), n1);
  r2 = min(ceil(sqrt((n2 / n1) * rsamp)), n2);
  r1 = max(1, r1);
  r2 = max(1, r2);
  G1 = randn(n1, r1);
  G2 = randn(n2, r2);
  PsiG = Psi * G1;
  PhiG = Phi * G2;
  BG = local_prod_states(PsiG, PhiG);
  [~, ~, e] = qr(BG', 0);
  clear BG;
  if numel(e) < rk
    error('isdf:gen_coeff_qrcp:Pivot', ...
      'QRCP returned only %d pivots; need rk=%d.', numel(e), rk);
  end
  ind_mu = int32(e(1:rk));
  ind_mu = ind_mu(:);

  output.msg('rs', 'gen_coeff_qrcp: id=%d  rk=%d  sketch=(%d,%d)  seed=%g', ...
    int32(id), rk, r1, r2, seed);

  % --- Fill sampling / coeff / bundle (fine-grid lin indices) ---
  Nmu = int32(numel(ind_mu));
  R_sampling_RLU = double(fft_data.Rgrid_RLU(ind_mu, :));

  nb = wf_data.nb;
  nspin = wf_data.nspin;
  coeff = zeros(double(Nmu), nb, nbz, nspin);
  for ikbz = 1:nbz
    ikibz = k_data.bz2ibz(ikbz, 1);
    ikrot = k_data.bz2rot(ikbz, 1);
    for ib = 1:nb
      wf = wave_functions.WF_apply_symm([ib, ikibz, ikrot, ispin]);
      coeff(:, ib, ikbz, ispin) = wf(ind_mu);
    end
  end

  % Local R_rot among selected points (identity fallback if orbit leaves the set).
  lin2loc = zeros(double(fft_data.nr), 1, 'int32');
  lin2loc(double(ind_mu)) = int32(1:double(Nmu));
  R_rot_coarse = zeros(double(Nmu), symm_data.nsym, 'int32');
  n_open = 0;
  if isempty(fft_data.R_rot) || size(fft_data.R_rot, 2) < symm_data.nsym
    R_rot_coarse(:, :) = repmat(int32(1:double(Nmu)).', 1, symm_data.nsym);
    n_open = double(Nmu) * double(symm_data.nsym);
  else
    for isym = 1:symm_data.nsym
      rot_lin = double(fft_data.R_rot(double(ind_mu), isym));
      loc = lin2loc(rot_lin);
      miss = loc <= 0;
      n_open = n_open + nnz(miss);
      loc(miss) = int32(find(miss));  % identity on rows that leave the set
      R_rot_coarse(:, isym) = loc;
    end
  end
  if n_open > 0
    output.warn( ...
      'gen_coeff_qrcp: QRCP set is not symmetry-closed (%d mappings left the set; used identity). Prefer adaptive refine.', ...
      n_open);
  end

  isdf_data.nisdf = Nmu;
  isdf_data.N_coarse = Nmu;
  isdf_data.N_extra = int32(0);
  isdf_data.coeff_seper = coeff;
  isdf_data.fftgrid_c = int32(fft_data.fftgrid(:).');
  isdf_data.R_sampling_RLU = R_sampling_RLU;
  isdf_data.interp_scheme = "qrcp";
  isdf_data.R_rot_coarse = R_rot_coarse;
  isdf_data.R_rot_extra = int32(zeros(0, 0));

  bundle_struct = struct();
  bundle_struct.N_bundle = isdf_data.N_coarse;
  bundle_struct.N_coarse = isdf_data.N_coarse;
  bundle_struct.N_sampling = isdf_data.N_coarse;
  bundle_struct.R_grid_bundle = R_sampling_RLU;
  bundle_struct.R_rot_in_bundle = R_rot_coarse;
  bundle_struct.WF_bundle = coeff;
  bundle_struct.sampling2bundle = int32((1:double(isdf_data.N_coarse)).');
  bundle_struct.fine_grid_lin = ind_mu;
  isdf_data.bundle_struct = bundle_struct;

  isdf.save2mod(isdf_data, id);
end

function P = local_prod_states(A, B)
% Element-wise products of columns: P(:, i+(j-1)*nA) = A(:,i) .* B(:,j).
  [m1, nA] = size(A);
  [m2, nB] = size(B);
  if m1 ~= m2
    error('isdf:gen_coeff_qrcp:prod_states', 'Row dimensions of A and B do not match.');
  end
  P = zeros(m1, nA * nB);
  for j = 1:nB
    for i = 1:nA
      P(:, i + (j - 1) * nA) = A(:, i) .* B(:, j);
    end
  end
end
