% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/19

% ISDF_COEFF.ISDF_APPLY_SYMM_ON_COARSE 鈥?package call: isdftest.coeff.isdf_apply_symm_on_coarse(...)
function cout = isdf_apply_symm_on_coarse(id, c_in, isym)
% Apply symmetry isym to values on coarse grid points (permute + optional conj for t-rev sector).
  persistent firsttime N_MAX nsym is_t_rev inv_rot_index R_rot

  if nargin == 1 && (ischar(id) || (isstring(id) && isscalar(id)))
    cmd = lower(string(id));
    if cmd == "reset"
      firsttime = false(N_MAX, 1);
      N_MAX = [];
      nsym = [];
      is_t_rev = [];
      inv_rot_index = [];
      R_rot = [];
      cout = [];
      return;
    end
  end

  if isempty(N_MAX)
    N_MAX = isdftest.isdftest_nmax();
    firsttime = true(N_MAX, 1);
    symm_data = symmetry.get();
    % We are facing the same system, so the symmetry properties are the same.
    nsym = symm_data.nsym;
    is_t_rev = symm_data.is_t_rev;
    inv_rot_index = symm_data.inv_rot_index;
    R_rot = cell(N_MAX, 1);
  end

  if (firsttime(id))
    firsttime(id) = false;
    isdf_data = isdftest.get(id);
    R_rot{id} = isdf_data.R_rot_coarse;
  end

  %
  isconj = false;
  if isym > nsym / (1 + is_t_rev)
    isym = isym - nsym / (1 + is_t_rev);
    isconj = true;
  end
  inv_isym = inv_rot_index(isym);
  
  
  ind = R_rot{id}(:, inv_isym);
  cout = c_in(ind, :);
  %
  if isconj 
    cout = conj(cout);
  end
  %
end
