% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function R_rot_coarse = isdf_coarse_R_rot(fftgrid_c)
% Symmetry index map on a regular coarse FFT RLU lattice (nr x nsym int32).

  fftgrid_c = int32(fftgrid_c(:).');
  if numel(fftgrid_c) ~= 3 || any(fftgrid_c < 1)
    error('isdf:isdf_coarse_R_rot:BadGrid', 'fftgrid_c must be 1x3 int32 with positive entries.');
  end
  gc_s = double(fftgrid_c);
  nr = prod(fftgrid_c);

  symm_data = symmetry.get();
  nsym = symm_data.nsym;
  is_t_rev = symm_data.is_t_rev;

  xq = 0:gc_s(1) - 1;
  yq = 0:gc_s(2) - 1;
  zq = 0:gc_s(3) - 1;
  [Xq, Yq, Zq] = ndgrid(xq, yq, zq);
  R_coarse_RLU = double([Xq(:), Yq(:), Zq(:)]);

  R_rot_coarse = zeros(nr, nsym, 'int32');
  % R_rot_coarse_wf = zeros(nr, nsym, 'int32');
  % R_rot_coarse_apply = zeros(nr, nsym, 'int32');

  for is = 1:nsym
    % if is <= nsym / (1 + is_t_rev)
    %   mtrx_RLU_R = symm_data.rot_mtrx_RLU_R(:, :, is);
    %   sign_factor = 1;
    % elseif is <= nsym
    %   mtrx_RLU_R = symm_data.rot_mtrx_RLU_R(:, :, is);
    %   sign_factor = -1;
    % end
    % M2 = double(sign_factor * mtrx_RLU_R);
    M2 = symm_data.rot_mtrx_RLU_R(:, :, is);

    M2_r_RLU = (R_coarse_RLU * M2);
    if norm(M2_r_RLU - round(M2_r_RLU)) > double(1e-4)
      error('isdf:isdf_coarse_R_rot:NonIntegerMap', ...
        'Non-integer mapping under rotation; check rot_mtrx_RLU_R and fftgrid_c.');
    end
    M2_r_RLU = round(M2_r_RLU);
    iv_mod = int32(mod(M2_r_RLU + gc_s, gc_s));
    dm = double(iv_mod);
    g1 = gc_s(1);
    g12 = gc_s(1) * gc_s(2);
    i4 = 1 + dm(:, 1) + dm(:, 2) * g1 + dm(:, 3) * g12;
    R_rot_coarse(:, is) = int32(round(i4));

    % if is > nsym / (1 + is_t_rev)
    %   M2 = double(mtrx_RLU_R);
    %   M2_r_RLU = (R_coarse_RLU * M2);
    %   iv_mod = int32(mod(M2_r_RLU + gc_s, gc_s));
    %   dm = double(iv_mod);
    %   g1 = gc_s(1);
    %   g12 = gc_s(1) * gc_s(2);
    %   i4 = 1 + dm(:, 1) + dm(:, 2) * g1 + dm(:, 3) * g12;
    % end
    % R_rot_coarse_apply(:, is) = int32(round(i4));
  end
end
