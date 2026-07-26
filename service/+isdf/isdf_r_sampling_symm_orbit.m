function [N_r_orbit, Nr, Rgrid_RLU, Rgrid_rep_RLU, ir2rep, ir2rot] = isdf_r_sampling_symm_orbit(id)
%ISDF_R_SAMPLING_SYMM_ORBIT Build symmetry-orbit mapping on ISDF sampling points.
%
% Usage:
%   isdf_r_sampling_symm_orbit()      % use current ISDF id in manager
%   isdf_r_sampling_symm_orbit(id)    % use explicit ISDF id
%
% Output convention:
%   N_r_orbit : number of orbits in ISDF R-sampling grid
%   Nr        : number of points in R_sampling_RLU
%   Rgrid_RLU : same as isdf_data.R_sampling_RLU
%   Rgrid_rep_RLU : representative points, size [N_r_orbit, 3]
%   ir2rep    : representative index for each sampling point
%   ir2rot    : rotation index for each sampling point
%
% Reconstruction relation:
%   Rgrid_RLU(ir, :) = Rgrid_rep_RLU(ir2rep(ir), :) * rot_mtrx_RLU_R(:, :, ir2rot(ir))

  if nargin < 1
    isdf_data = isdf.get();
  else
    isdf_data = isdf.get(id);
  end

  fft_data = FFT.get();
  symm_data = symmetry.get();

  Rgrid_RLU = double(isdf_data.R_sampling_RLU);
  Nr = int32(size(Rgrid_RLU, 1));
  fftgrid = int32(fft_data.fftgrid);
  fftgrid_d = double(fftgrid);
  rot_mtrx_RLU_R = double(symm_data.rot_mtrx_RLU_R);
  nsym = int32(size(rot_mtrx_RLU_R, 3));
  is_t_rev = symm_data.is_t_rev;

  if Nr <= 0
    error('isdf_r_sampling_symm_orbit:EmptySampling', 'isdf_data.R_sampling_RLU is empty.');
  end
  if nsym <= 0
    error('isdf_r_sampling_symm_orbit:NoSymmetry', 'symmetry.rot_mtrx_RLU_R is empty.');
  end

  tol = 1e-4;
  R_mod_all = local_mod_rows(double(Rgrid_RLU), fftgrid_d);

  % Build point-index map for O(1) lookup on rotated images.
  point2ir = containers.Map('KeyType', 'char', 'ValueType', 'any');
  for ir = 1:Nr
    key = local_key(R_mod_all(ir, :), tol);
    if isKey(point2ir, key)
      point2ir(key) = [point2ir(key), int32(ir)];
    else
      point2ir(key) = int32(ir);
    end
  end

  ir2rep = int32(zeros(Nr, 1));
  ir2rot = int32(zeros(Nr, 1));
  Rgrid_rep_RLU = double(zeros(0, 3));
  rep_count = int32(0);

  for ir = 1:Nr
    if ir2rep(ir) ~= 0
      continue;
    end

    rep_count = rep_count + 1;
    r_rep = double(Rgrid_RLU(ir, :));
    Rgrid_rep_RLU(rep_count, :) = double(r_rep); %#ok<AGROW>

    for irot = 1 : nsym / (1 + is_t_rev)
      r_img = local_mod_rows(r_rep * rot_mtrx_RLU_R(:, :, irot), fftgrid_d);
      key_img = local_key(r_img, tol);
      if ~isKey(point2ir, key_img)
        continue;
      end

      ir_list = point2ir(key_img);
      if ~isa(ir_list, 'int32')
        ir_list = int32(ir_list);
      end

      for k = 1:numel(ir_list)
        ir_img = ir_list(k);
        if local_periodic_dist(r_img, R_mod_all(ir_img, :), fftgrid_d) > tol
          continue;
        end
        if ir2rep(ir_img) == 0
          ir2rep(ir_img) = rep_count;
          ir2rot(ir_img) = irot;
        end
      end
    end
  end

  if any(ir2rep == 0) || any(ir2rot == 0)
    error('isdf_r_sampling_symm_orbit:IncompleteCover', ...
      'R_sampling_RLU is not closed under symmetry; some points are uncovered.');
  end

  N_r_orbit = rep_count;

  % Consistency check against reconstruction relation.
  for ir = 1:Nr
    irep = ir2rep(ir);
    irot = ir2rot(ir);
    r_rhs = local_mod_rows(double(Rgrid_rep_RLU(irep, :)) * rot_mtrx_RLU_R(:, :, irot), fftgrid_d);
    r_lhs = R_mod_all(ir, :);
    if local_periodic_dist(r_rhs, r_lhs, fftgrid_d) > tol
      error('isdf_r_sampling_symm_orbit:RelationMismatch', ...
        'Orbit relation mismatch at ir=%d.', ir);
    end
  end
end

function rmod = local_mod_rows(r, period)
  rmod = mod(r, period);
end

function d = local_periodic_dist(a, b, period)
  delta = abs(a - b);
  delta = min(delta, abs(period - delta));
  d = max(delta);
end

function key = local_key(rmod, tol)
  q = int64(round(rmod ./ tol));
  key = sprintf('%d_%d_%d', q(1), q(2), q(3));
end
