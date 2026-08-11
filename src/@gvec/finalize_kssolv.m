% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2025/07/18 ZZ

function gvec_out = finalize_kssolv(mol, ecut, gvec)
  grid = Ggrid(mol, ecut);
  ng = grid.ng;
  coulG = [grid.gkx, grid.gky, grid.gkz] * mol.supercell / (2*pi);
  coulG = round(coulG);
  gvec.components = coulG;
  % gvec.index_vec = zeros(mol.n1 * mol.n2 * mol.n3, 1);
  % gvec.index_vec(grid.idxnz) = (1:ng)'
  gvec.ng = ng;
  gvec.nfftgridpts = prod([mol.n1, mol.n2, mol.n3]);
  gvec.fftgrid = [mol.n1, mol.n2, mol.n3];
  return; 
end
