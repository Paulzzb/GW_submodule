% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/11

function varargout = SC_ISDF(source, ratio, varargin)
%SC_ISDF  Replicate unit-cell ISDF sampling grid to a k-supercell.
%
%   Suppose the supercell is a k = [k1, k2, k3] replication of the unit cell
%   with a consistent real-space FFT grid: if the unit cell uses [N1, N2, N3],
%   the supercell uses [N1*k1, N2*k2, N3*k3].
%
%   id_out = isdf.SC_ISDF(id, ratio)
%   id_out = isdf.SC_ISDF(mat_path, ratio)
%   [id_out, info] = isdf.SC_ISDF(..., 'SavePath', fpath)
%   [id_out, info] = isdf.SC_ISDF(..., 'UcFftgrid', [N1 N2 N3])
%
%   Steps (README):
%     a. Read isdf*.mat (pool id or .mat containing isdf_data)
%     b. Read the grid after adaptiveisdf (R_sampling_RLU)
%     c. Optional QwQ: if (nisdf_uc * prod(ratio)) / sqrt(Nn1*Nn2) exceeds
%        TargetIsdfRatio, randomly subsample UC points before replication
%     d. Duplicate the grid to each unit cell inside the supercell
%     e. Recalculate linear indices on the supercell FFT grid and save
%
%   Inputs:
%     source  - ISDF pool id (numeric) or path to .mat with variable isdf_data
%     ratio   - 1x3 positive integers [k1, k2, k3]
%
%   Name-value options:
%     SavePath   - if set, write isdf_data (and metadata) to this .mat path
%     UcFftgrid        - unit-cell FFT grid [N1 N2 N3]; default: FFT.get().fftgrid ./ ratio
%     TargetIsdfRatio  - max nisdf/sqrt(Nn1*Nn2) before replication (0 = off)
%     SystemCfg        - config.SYSTEM for nrange (required when TargetIsdfRatio > 0)
%     IdOut            - target pool id; default: allocate a new slot
%     SourceId         - id_coarse metadata when loading from .mat (optional)
%
%   Outputs:
%     id_out - pool id of the supercell ISDF object
%     info   - struct with uc_fftgrid, sc_fftgrid, ratio, nisdf_uc, nisdf_sc, save_path

  p = inputParser;
  p.addRequired('source');
  p.addRequired('ratio');
  p.addParameter('SavePath', '', @(x) ischar(x) || isstring(x));
  p.addParameter('UcFftgrid', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
  p.addParameter('TargetIsdfRatio', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
  p.addParameter('SystemCfg', [], @(x) isempty(x) || isstruct(x));
  p.addParameter('IdOut', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
  p.addParameter('SourceId', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
  p.parse(source, ratio, varargin{:});
  opts = p.Results;

  ratio = int32(round(double(ratio(:)).'));
  if numel(ratio) ~= 3 || any(ratio < 1)
    error('isdf:SC_ISDF:BadRatio', 'ratio must be three positive integers [k1, k2, k3].');
  end

  [isdf_data, source_meta] = sc_isdf_load_source(source);
  uc_fftgrid = sc_isdf_resolve_uc_fftgrid(opts.UcFftgrid, ratio);
  sc_fftgrid = int32(double(uc_fftgrid) .* double(ratio));

  fft_data = FFT.get();
  if ~isequal(int32(fft_data.fftgrid(:).'), sc_fftgrid(:).')
    error('isdf:SC_ISDF:FftgridMismatch', ...
      ['FFT.get().fftgrid = [%s] but supercell grid requires [%s]. ', ...
       'Run FFT.driver on the supercell first or pass UcFftgrid explicitly.'], ...
      num2str(double(fft_data.fftgrid)), num2str(double(sc_fftgrid)));
  end

  Nmu_uc = double(isdf_data.nisdf);
  if Nmu_uc < 1
    error('isdf:SC_ISDF:EmptyIsdf', 'Input ISDF has nisdf < 1.');
  end
  if isempty(isdf_data.R_sampling_RLU)
    error('isdf:SC_ISDF:EmptyGrid', 'Input ISDF R_sampling_RLU is empty.');
  end

  N_coarse_uc = double(isdf_data.N_coarse);
  N_extra_uc = double(isdf_data.N_extra);
  if N_coarse_uc <= 0 && N_extra_uc <= 0
    N_coarse_uc = Nmu_uc;
    N_extra_uc = 0;
  elseif N_coarse_uc + N_extra_uc ~= Nmu_uc
    error('isdf:SC_ISDF:NmuSplit', ...
      'N_coarse (%d) + N_extra (%d) != nisdf (%d).', ...
      int32(N_coarse_uc), int32(N_extra_uc), int32(isdf_data.nisdf));
  end

  Rs_uc = double(isdf_data.R_sampling_RLU(1:Nmu_uc, :));
  qwq_info = struct('applied', false);
  target_ratio = double(opts.TargetIsdfRatio);
  if target_ratio > 0
    if isempty(opts.SystemCfg)
      error('isdf:SC_ISDF:QwQConfig', ...
        'SystemCfg (config.SYSTEM) is required when TargetIsdfRatio > 0.');
    end
    [Rs_uc, N_coarse_uc, N_extra_uc, qwq_info] = sc_isdf_qwq_trim( ...
      Rs_uc, N_coarse_uc, N_extra_uc, target_ratio, ...
      char(string(isdf_data.desc)), opts.SystemCfg, ratio);
    Nmu_uc = size(Rs_uc, 1);
  end
  [R_sc, lin_sc, N_coarse_sc, N_extra_sc] = sc_isdf_duplicate_sampling( ...
    Rs_uc, N_coarse_uc, N_extra_uc, uc_fftgrid, ratio);

  isdf_sc = sc_isdf_build_supercell_object(isdf_data, R_sc, lin_sc, ...
    uc_fftgrid, sc_fftgrid, ratio, N_coarse_sc, N_extra_sc, fft_data);

  id_out = sc_isdf_write_pool(isdf_sc, opts.IdOut);
  isdf_sc.id = id_out;

  save_path = char(string(opts.SavePath));
  if strlength(string(save_path)) > 0
    sc_isdf_save_mat(save_path, isdf_sc, source_meta, ratio, uc_fftgrid, sc_fftgrid, id_out);
  else
    save_path = '';
  end

  info = struct();
  info.id_out = id_out;
  info.ratio = double(ratio);
  info.uc_fftgrid = double(uc_fftgrid);
  info.sc_fftgrid = double(sc_fftgrid);
  info.nisdf_uc = Nmu_uc;
  if qwq_info.applied
    info.nisdf_uc_before = qwq_info.nisdf_before;
  else
    info.nisdf_uc_before = Nmu_uc;
  end
  info.nisdf_sc = double(isdf_sc.nisdf);
  info.qwq = qwq_info;
  info.N_coarse_uc = N_coarse_uc;
  info.N_extra_uc = N_extra_uc;
  info.N_coarse_sc = N_coarse_sc;
  info.N_extra_sc = N_extra_sc;
  info.save_path = save_path;
  info.source = source_meta;

  fprintf(['SC_ISDF: %s -> pool id=%d, nisdf %d -> %d, ', ...
    'uc fft [%s], sc fft [%s], ratio [%s]\n'], ...
    char(string(source_meta.label)), double(id_out), Nmu_uc, double(isdf_sc.nisdf), ...
    num2str(double(uc_fftgrid)), num2str(double(sc_fftgrid)), num2str(double(ratio)));

  if nargout >= 1
    varargout{1} = id_out;
  end
  if nargout >= 2
    varargout{2} = info;
  end
end

function [isdf_data, meta] = sc_isdf_load_source(source)
  meta = struct('kind', '', 'label', '', 'path', '', 'id', []);

  if isnumeric(source) && isscalar(source)
    id = int32(source);
    isdf_data = isdf.get(id);
    meta.kind = 'pool';
    meta.label = sprintf('pool id=%d', double(id));
    meta.id = id;
    return;
  end

  fpath = char(string(source));
  if exist(fpath, 'file') ~= 2
    error('isdf:SC_ISDF:BadSource', ...
      'source must be a pool id or an existing .mat path (got ''%s'').', fpath);
  end
  S = load(fpath);
  if ~isfield(S, 'isdf_data')
    error('isdf:SC_ISDF:BadMat', 'File %s must contain variable isdf_data.', fpath);
  end
  isdf_data = S.isdf_data;
  if ~isa(isdf_data, 'isdf.base.isdf_m')
    error('isdf:SC_ISDF:BadMatType', ...
      'isdf_data in %s must be isdf.base.isdf_m (got %s).', fpath, class(isdf_data));
  end
  meta.kind = 'mat';
  meta.label = fpath;
  meta.path = fpath;
  if isfield(S, 'id_adaptive')
    meta.id_adaptive = S.id_adaptive;
  end
  if isfield(S, 'id_coarse')
    meta.id_coarse = S.id_coarse;
  end
  if isfield(S, 'isdf_desc')
    meta.isdf_desc = char(string(S.isdf_desc));
  end
end

function uc_fftgrid = sc_isdf_resolve_uc_fftgrid(uc_override, ratio)
  if ~isempty(uc_override)
    uc_fftgrid = int32(round(double(uc_override(:)).'));
    if numel(uc_fftgrid) ~= 3 || any(uc_fftgrid < 1)
      error('isdf:SC_ISDF:BadUcFftgrid', 'UcFftgrid must be three positive integers.');
    end
    return;
  end

  fft_data = FFT.get();
  sc_fftgrid = int32(fft_data.fftgrid(:).');
  if numel(sc_fftgrid) ~= 3
    error('isdf:SC_ISDF:BadFftgrid', 'FFT.get().fftgrid must have 3 components.');
  end
  uc_fftgrid = int32(round(double(sc_fftgrid) ./ double(ratio)));
  if ~isequal(int32(double(uc_fftgrid) .* double(ratio)), sc_fftgrid)
    error('isdf:SC_ISDF:NotCommensurate', ...
      'FFT grid [%s] is not commensurate with ratio [%s]. Pass UcFftgrid explicitly.', ...
      num2str(double(sc_fftgrid)), num2str(double(ratio)));
  end
end

function [Rs_uc, N_coarse_uc, N_extra_uc, info] = sc_isdf_qwq_trim( ...
    Rs_uc, N_coarse_uc, N_extra_uc, target_ratio, desc, system_cfg, sc_ratio)

  info = struct('applied', false, 'target_ratio', target_ratio, ...
    'nisdf_before', size(Rs_uc, 1), 'nisdf_after', size(Rs_uc, 1), ...
    'nmu_target', NaN, 'ratio_before', NaN, 'ratio_after', NaN);

  [~, ~, Nn1, Nn2] = isdf.resolve_nrange(desc, system_cfg);
  denom = sqrt(double(Nn1) * double(Nn2));
  nrep = prod(double(sc_ratio(:)));
  % Cap UC points so post-replication SC ratio <= target_ratio:
  %   (nisdf_uc * nrep) / denom <= target_ratio
  nmu_target = max(1, floor(target_ratio * denom / nrep));
  Nmu = size(Rs_uc, 1);
  ratio_before = (Nmu * nrep) / denom;

  info.nmu_target = nmu_target;
  info.ratio_before = ratio_before;
  info.ratio_after = ratio_before;

  if Nmu <= nmu_target
    return;
  end

  keep = sort(randperm(Nmu, nmu_target));
  Rs_uc = Rs_uc(keep, :);
  N_coarse_uc = sum(keep <= N_coarse_uc);
  N_extra_uc = nmu_target - N_coarse_uc;

  info.applied = true;
  info.nisdf_before = Nmu;
  info.nisdf_after = nmu_target;
  info.ratio_after = (nmu_target * nrep) / denom;

  fprintf(['SC_ISDF QwQ: desc=%s nisdf %d -> %d ', ...
    '(SC ISDF ratio %.3f -> %.3f, cap %.3f)\n'], ...
    desc, Nmu, nmu_target, ratio_before, info.ratio_after, target_ratio);
end

function [R_sc, lin_sc, N_coarse_sc, N_extra_sc] = sc_isdf_duplicate_sampling( ...
    Rs_uc, N_coarse_uc, N_extra_uc, uc_fftgrid, ratio)

  uc = double(uc_fftgrid(:)).';
  k = double(ratio(:)).';
  sc = uc .* k;
  nrep = prod(k);

  blocks = {Rs_uc(1:N_coarse_uc, :), Rs_uc(N_coarse_uc + (1:N_extra_uc), :)};
  block_sizes = [N_coarse_uc, N_extra_uc];
  R_parts = cell(1, 2);
  lin_parts = cell(1, 2);

  for ib = 1:2
    if block_sizes(ib) < 1
      R_parts{ib} = zeros(0, 3);
      lin_parts{ib} = zeros(0, 1);
      continue;
    end
    Rb = blocks{ib};
    R_out = zeros(block_sizes(ib) * nrep, 3);
    lin_out = zeros(block_sizes(ib) * nrep, 1);
    row = 0;
    for i3 = 0:k(3) - 1
      for i2 = 0:k(2) - 1
        for i1 = 0:k(1) - 1
          offset = [i1 * uc(1), i2 * uc(2), i3 * uc(3)];
          R_shift = mod(round(Rb + offset), sc);
          idx = row + (1:block_sizes(ib));
          R_out(idx, :) = R_shift;
          lin_out(idx) = sc_isdf_rlu_to_lin(R_shift, sc);
          row = row + block_sizes(ib);
        end
      end
    end
    R_parts{ib} = R_out;
    lin_parts{ib} = lin_out;
  end

  R_sc = [R_parts{1}; R_parts{2}];
  lin_sc = [lin_parts{1}; lin_parts{2}];
  N_coarse_sc = block_sizes(1) * nrep;
  N_extra_sc = block_sizes(2) * nrep;
end

function lin = sc_isdf_rlu_to_lin(R_rlu, fftgrid)
  g = double(fftgrid(:)).';
  iv = mod(round(double(R_rlu)), g);
  lin = 1 + iv(:, 1) + iv(:, 2) * g(1) + iv(:, 3) * g(1) * g(2);
  lin = int32(lin);
end

function isdf_sc = sc_isdf_build_supercell_object(isdf_uc, R_sc, lin_sc, ...
    uc_fftgrid, sc_fftgrid, ratio, N_coarse_sc, N_extra_sc, fft_data)

  isdf_sc = isdf_uc;
  nmu_sc = size(R_sc, 1);
  isdf_sc.nisdf = int32(nmu_sc);
  isdf_sc.R_sampling_RLU = R_sc;
  isdf_sc.N_coarse = int32(N_coarse_sc);
  isdf_sc.N_extra = int32(N_extra_sc);

  if ~isempty(isdf_uc.fftgrid_c)
    isdf_sc.fftgrid_c = int32(round(double(isdf_uc.fftgrid_c(:)).' .* double(ratio)));
  else
    isdf_sc.fftgrid_c = int32(zeros(1, 0));
  end

  isdf_sc.R_rot_extra = sc_isdf_build_R_rot_extra(lin_sc, fft_data);
  if double(isdf_sc.N_coarse) > 0
    isdf_sc.R_rot_coarse = isdf_sc.R_rot_extra(1:double(isdf_sc.N_coarse), :);
  else
    isdf_sc.R_rot_coarse = int32(zeros(0, 0));
  end

  isdf_sc.coeff_seper = zeros(0, 0, 0, 0);
  isdf_sc.tildeVq = zeros(0, 0, 0, 0);
  isdf_sc.helperqG = zeros(0, 0, 0, 0);
  isdf_sc.CCHq = zeros(0, 0, 0);
  isdf_sc.CCHq_trunc_factors = {};
  isdf_sc.tmp = [];

  desc_uc = char(string(isdf_uc.desc));
  isdf_sc.desc = string(sprintf('%s_sc_%d_%d_%d', desc_uc, ratio(1), ratio(2), ratio(3)));
  if isdf_sc.interp_scheme == "coarse"
    isdf_sc.interp_scheme = "coarse_sc";
  elseif isdf_sc.interp_scheme == "adaptive"
    isdf_sc.interp_scheme = "adaptive_sc";
  else
    isdf_sc.interp_scheme = "supercell";
  end

  bs = struct();
  bs.N_sampling = int32(nmu_sc);
  bs.N_coarse = isdf_sc.N_coarse;
  bs.N_bundle = int32(0);
  bs.fine_grid_lin = int32(lin_sc(:));
  bs.sampling2bundle = int32((1:nmu_sc).');
  bs.R_grid_bundle = R_sc;
  bs.R_rot_in_bundle = int32(zeros(0, 0));
  bs.WF_bundle = zeros(0, 0, 0, 0);
  bs.sc_ratio = double(ratio);
  bs.uc_fftgrid = double(uc_fftgrid);
  bs.sc_fftgrid = double(sc_fftgrid);
  isdf_sc.bundle_struct = bs;
  isdf_sc.assigned = true;
end

function R_rot_extra = sc_isdf_build_R_rot_extra(lin_sc, fft_data)
  lin_sc = int32(lin_sc(:));
  Ntot = numel(lin_sc);
  if Ntot < 1
    R_rot_extra = int32(zeros(0, 0));
    return;
  end

  symm_data = symmetry.get();
  nsym = double(symm_data.nsym);
  fftgrid = double(fft_data.fftgrid(:)).';
  Rgrid = double(fft_data.Rgrid_RLU);

  R_rot_extra = zeros(Ntot, nsym, 'int32');
  for is = 1:nsym
    mtrx_RLU_R = symm_data.rot_mtrx_RLU_R(:, :, is);
    M2_r_RLU = Rgrid(lin_sc, :) * mtrx_RLU_R;
    if max(abs(M2_r_RLU - round(M2_r_RLU)), [], 'all') > 1e-3
      error('isdf:SC_ISDF:NonIntegerRot', ...
        'Non-integer rotation map on supercell FFT grid (symmetry %d).', is);
    end
    iv_mod = int32(mod(round(M2_r_RLU) + fftgrid, fftgrid));
    i4 = 1 + iv_mod(:, 1) + iv_mod(:, 2) * fftgrid(1) ...
      + iv_mod(:, 3) * fftgrid(1) * fftgrid(2);
    R_rot_extra(:, is) = int32(i4);
  end
end

function id_out = sc_isdf_write_pool(isdf_sc, id_out_opt)
  if ~isempty(id_out_opt)
    id_out = int32(id_out_opt);
  else
    id_out = isdf.isdf_add(isdf_sc.desc);
  end
  isdf_sc.id = id_out;
  isdf.save2mod(isdf_sc, id_out);
end

function sc_isdf_save_mat(fpath, isdf_data, source_meta, ratio, uc_fftgrid, sc_fftgrid, id_out)
  parentDir = fileparts(fpath);
  if ~isempty(parentDir) && exist(parentDir, 'dir') ~= 7
    mkdir(parentDir);
  end
  isdf_desc = char(string(isdf_data.desc));
  sc_ratio = double(ratio);
  uc_fftgrid = double(uc_fftgrid);
  sc_fftgrid = double(sc_fftgrid);
  id_sc = id_out;
  source = source_meta;
  save(fpath, 'isdf_data', 'isdf_desc', 'sc_ratio', 'uc_fftgrid', 'sc_fftgrid', 'id_sc', 'source', '-v7.3');
  fprintf('SC_ISDF: wrote %s\n', fpath);
end
