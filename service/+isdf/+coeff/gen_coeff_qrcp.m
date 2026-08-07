% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ
%
% QRCP index route ported from kssolvGW .../ISDF/isdf_indices.m (case 'qrcp').

function idx_mu = gen_coeff_qrcp(cfg, id)
%GEN_COEFF_QRCP  Select ISDF sampling points by randomized QRCP.
%
%   idx_mu = isdf.coeff.gen_coeff_qrcp(cfg, id)
%
% Returns fine-grid linear indices only. Caller fills the ISDF slot via
% isdf.rsymm.init_from_indices(id, idx_mu).
%
% Algorithm (legacy isdf_indices / exxmethod='qrcp'):
%   1. Build real-space Psi/Phi on the fine FFT grid from nrange1/nrange2 (all BZ k).
%   2. Rank rk = isdf_data.nisdf (set by driver from isdf_ratio * sqrt(N1*N2)).
%   3. Gaussian-sketch Phi/Psi, form product states, QR with column pivoting;
%      take the first rk pivoted row indices as fine-grid sampling points.
%
% cfg may supply:
%   .seed  - RNG seed (default 0, same as config.ISDF.seed default)

  if nargin < 1 || isempty(cfg)
    cfg = struct();
  end
  if nargin < 2 || isempty(id)
    output.err('ISDF id is required.');
  end

  fft_data = FFT.get();
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  isdf_data = isdf.get(id);

  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    output.err( ...
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
    output.err('Psi/Phi row dimensions do not match.');
  end

  rk = double(isdf_data.nisdf);
  if ~(isfinite(rk) && rk >= 1)
    output.err( ...
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
  BG = prod_states(PsiG, PhiG);
  [~, ~, e] = qr(BG', 0);
  clear BG;
  if numel(e) < rk
    output.err( ...
      'QRCP returned only %d pivots; need rk=%d.', numel(e), rk);
  end
  idx_mu = int32(e(1:rk));
  idx_mu = idx_mu(:);

  output.msg('rs', 'gen_coeff_qrcp: id=%d  rk=%d  sketch=(%d,%d)  seed=%g', ...
    int32(id), rk, r1, r2, seed);
end

function P = prod_states(A, B)
% Element-wise products of columns: P(:, i+(j-1)*nA) = A(:,i) .* B(:,j).
  [m1, nA] = size(A);
  [m2, nB] = size(B);
  if m1 ~= m2
    output.err('Row dimensions of A and B do not match.');
  end
  P = zeros(m1, nA * nB);
  for j = 1:nB
    for i = 1:nA
      P(:, i + (j - 1) * nA) = A(:, i) .* B(:, j);
    end
  end
end
