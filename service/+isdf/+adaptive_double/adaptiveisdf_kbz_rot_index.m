% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/16

function ikrotbz = adaptiveisdf_kbz_rot_index(ikbz, irot)
% adaptiveisdf_kbz_rot_index  Map full-BZ k index under a crystal symmetry in RLU.
%
%   ikrotbz = adaptiveisdf_kbz_rot_index(ikbz, irot)
%
% Inputs:
%   ikbz  in 1 .. Nbz   (full BZ index, same ordering as lattice.k / kptbz_RLU rows)
%   irot  in 1 .. nsym (symmetry.manager / symm_m index, including time-reversal sector)
%
% Output:
%   ikrotbz  such that (row vectors, reciprocal lattice / RLU coordinates)
%
%     kptbz_RLU(ikrotbz,:)  ==  kptbz_RLU(ikbz,:) * rot_mtrx_RLU_G(:,:,irot)
%
%   up to a reciprocal-lattice vector (components differ by integers), with
%   tolerance r_lat.tol (same spirit as KPT_expand).
%
% Requires lattice.manager('k','get') and symmetry.manager('get') initialized.

  k_data = lattice.manager('k', 'get');
  symm_data = symmetry.manager('get');
  r_lat = lattice.manager('r_lat', 'get');

  nbz = double(k_data.nbz);
  nsym = double(symm_data.nsym);
  validateattributes(ikbz, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'ikbz');
  validateattributes(irot, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'irot');
  if ikbz > nbz
    error('adaptiveisdf:kbz_rot_index:ikbz', 'ikbz=%d exceeds Nbz=%d.', ikbz, nbz);
  end
  if irot > nsym
    error('adaptiveisdf:kbz_rot_index:irot', 'irot=%d exceeds nsym=%d.', irot, nsym);
  end

  k_row = double(k_data.kptbz_RLU(ikbz, :));
  R = double(symm_data.rot_mtrx_RLU_G(:, :, irot));
  k_rot = k_row * R;

  tol = double(r_lat.tol);
  if ~(tol > 0)
    tol = 1e-5;
  end

  for j = 1:nbz
    k_j = double(k_data.kptbz_RLU(j, :));
    d = k_rot - k_j;
    d = d - round(d);
    if norm(d) < tol
      ikrotbz = int32(j);
      return
    end
  end

  error('adaptiveisdf:kbz_rot_index:noMatch', ...
    'No full-BZ k-point matches k(%d,:)*R(:,:,%d) up to integer G* (Nbz=%d).', ikbz, irot, nbz);
end
