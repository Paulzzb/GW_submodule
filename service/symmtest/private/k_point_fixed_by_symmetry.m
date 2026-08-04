% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/19 ZZ

function tf = k_point_fixed_by_symmetry(ik, isym, kptbz, rot_mtrx_RLU_G, tol)
% True iff S*k = k in fractional reciprocal coordinates (mod integer G*).
S = double(rot_mtrx_RLU_G(:, :, isym));
k = double(kptbz(ik, :));
diff = k * S - k;
err = norm(diff - round(diff));
tf = err <= tol;
end
