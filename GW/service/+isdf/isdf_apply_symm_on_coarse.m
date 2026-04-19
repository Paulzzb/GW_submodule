% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function cout = isdf_apply_symm_on_coarse(id, c_in, isym)
% Apply symmetry isym to values on coarse grid points (permute + optional conj for t-rev sector).
  persistent firsttime nsym is_t_rev inv_rot_index R_rot_coarse

  if nargin == 1 && (ischar(id) || (isstring(id) && isscalar(id)))
    cmd = lower(string(id));
    if cmd == "reset"
      firsttime = [];
      nsym = [];
      is_t_rev = [];
      inv_rot_index = [];
      R_rot_coarse = [];
      cout = [];
      return;
    end
  end

  if isempty(firsttime)
    firsttime = true;
  end
  %
  if firsttime
    symm_data = symmetry.get();
    nsym = symm_data.nsym;
    is_t_rev = symm_data.is_t_rev;
    inv_rot_index = symm_data.inv_rot_index;
    %
    isdf_data = isdf.get(id);
    R_rot_coarse = isdf_data.R_rot_coarse;
    %
    firsttime = false;
  end
  %
  if isym == 1
    cout = c_in;
    return;
  end
  %
  inv_isym = inv_rot_index(isym);
  ind = R_rot_coarse(:, inv_isym);
  cout = c_in(ind, :);
  %
  if isym > nsym / (1 + is_t_rev)
    cout = conj(cout);
  end
  %
end
