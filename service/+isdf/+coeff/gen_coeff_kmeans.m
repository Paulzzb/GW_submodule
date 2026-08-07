% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ
%
% K-means index route ported from kssolvGW .../ISDF/isdf_indices.m (case 'kmeans').

function idx_mu = gen_coeff_kmeans(cfg, id)
%GEN_COEFF_KMEANS  Select ISDF sampling points by weighted k-means.
%
%   idx_mu = isdf.coeff.gen_coeff_kmeans(cfg, id)
%
% Returns fine-grid linear indices only. Caller fills the ISDF slot via
% isdf.rsymm.init_from_indices(id, idx_mu).
%
% Algorithm (legacy isdf_indices / exxmethod='kmeans'):
%   1. Build real-space Psi/Phi on the fine FFT grid from nrange1/nrange2 (all BZ k).
%   2. Rank rk = isdf_data.nisdf.
%   3. Build grid weight (add/prod/power), then k_means on weight.^2.
%
% cfg may supply (from config.ISDF):
%   .weight  - 'add' (default) | 'prod' | 'power'
%   .power   - used when weight='power'
%   .seed    - RNG seed (default 0)
%   .init    - k-means init (default 'random')

  if nargin < 1 || isempty(cfg)
    cfg = struct();
  end
  if nargin < 2 || isempty(id)
    output.err('ISDF id is required.');
  end

  fft_data = FFT.get();
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  d_lat = lattice.manager('d_lat', 'get');
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

  % --- Build Phi/Psi on the fine grid (same layout as gen_coeff_qrcp) ---
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

  if size(Psi, 1) ~= size(Phi, 1)
    output.err('Psi/Phi row dimensions do not match.');
  end

  rk = double(isdf_data.nisdf);
  if ~(isfinite(rk) && rk >= 1)
    output.err( ...
      'isdf_data.nisdf must be a positive rank before k-means (id=%d).', int32(id));
  end
  rk = min([rk, Nw]);
  rk = max(1, round(rk));

  weight_kind = 'add';
  if isstruct(cfg) && isfield(cfg, 'weight') && ~isempty(cfg.weight)
    weight_kind = lower(char(string(cfg.weight)));
  end
  switch weight_kind
    case 'prod'
      weight = sum(abs(Phi).^2, 2) .* sum(abs(Psi).^2, 2);
    case 'add'
      weight = sum(abs(Phi).^2, 2) + sum(abs(Psi).^2, 2);
    case 'power'
      if ~isfield(cfg, 'power') || isempty(cfg.power)
        output.err('weight=''power'' requires cfg.power.');
      end
      weight = (sum(abs(Phi).^2, 2) .* sum(abs(Psi).^2, 2)).^(double(cfg.power) / 2);
    case 'hf'
      output.err('weight=''hf'' is not supported in gen_coeff_kmeans yet.');
    otherwise
      output.err('Unknown weight kind ''%s''.', weight_kind);
  end
  weight_square = weight.^2;

  seed = 0;
  if isstruct(cfg) && isfield(cfg, 'seed') && ~isempty(cfg.seed)
    seed = double(cfg.seed);
  end
  init = 'random';
  if isstruct(cfg) && isfield(cfg, 'init') && ~isempty(cfg.init)
    init = char(string(cfg.init));
  end

  opt = struct();
  opt.seed = seed;
  opt.init = init;
  opt.fftgrid = double(fft_data.fftgrid(:).');
  opt.supercell = double(d_lat.a1a2a3);
  if ~isempty(d_lat.atom_pos)
    opt.atompos = reshape(double(d_lat.atom_pos), size(d_lat.atom_pos, 1), 3);
  else
    opt.atompos = [];
  end

  idx_mu = isdf.coeff.k_means(rk, weight_square, opt);
  idx_mu = int32(idx_mu(:));

  output.msg('rs', 'gen_coeff_kmeans: id=%d  rk=%d  weight=%s  init=%s  seed=%g', ...
    int32(id), rk, weight_kind, init, seed);
end
