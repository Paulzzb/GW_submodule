function X0 = genX0(mol,nXcols)
% GENX0 generates the initial wave functions.
%    X0 = GENX0(mol,ncols) generates the initial wave functions for
%    molecule mol with ncols columns.
%
%    BX0 = GENX0(cry,ncols) generates the initial Bloch wave functions for
%    crystal cry with ncols columns for each k-points.
%
%   See also Molecule, Crystal, Wavefun, BlochWavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nspin = mol.nspin;
lsda = mol.lsda;
noncolin = mol.noncolin;

if isa( mol, 'Crystal' )
    nkpts = mol.nkpts;
else
    nkpts = 1;
end

if nargin == 1
    if isa( mol, 'Crystal' )
        if noncolin
            nXcols = mol.nel*ones(nkpts,1);
        else
            nXcols = mol.nel/2*ones(nkpts*nspin,1);
        end
    else
        if noncolin
            nXcols = mol.nel;
        else
            nXcols = mol.nel/2;
        end
    end
end

n1   = mol.n1;
n2   = mol.n2;
n3   = mol.n3;

grid = Ggrid(mol);
idxnz = grid.idxnz;
npw = length(idxnz);
npol = noncolin + 1;

Qcell = cell(nkpts*(1+lsda),1);

psif = zeros(npw*npol,max(nXcols));
for ik = 1:nkpts
  for j = 1:nXcols(ik)
      for ipol = 1:npol
          psir = randn(n1,n2,n3);
          psig = reshape(fft3(psir),n1*n2*n3,1);
          idx = (1:npw) + (ipol-1)*npw;
          psif(idx,j) = psig(idxnz,1);
      end
  end
  % Cholesky is faster than QR, but is less stable.
  try
    Qcell{ik}=psif/chol(psif'*psif);
  catch
    [Qcell{ik},~]=qr(psif,0);
  end
end

if lsda
    for ik = 1:nkpts
        for j = 1:nXcols(ik)
            psir = randn(n1,n2,n3);
            psig = reshape(fft3(psir),n1*n2*n3,1);
            psif(:,j) = psig(idxnz,1);
        end
    end

    try
        Qcell{ik+nkpts}=psif/chol(psif'*psif);
    catch
        [Qcell{ik+nkpts},~]=qr(psif,0);
    end
end

if isa( mol, 'Crystal' )
    X0 = BlochWavefun(Qcell,n1,n2,n3,idxnz,mol.wks,nspin);
else
    if ~lsda
        X0 = Wavefun(Qcell{1},n1,n2,n3,idxnz);
    else
        X0 = cell(2,1);
        X0{1} = Wavefun(Qcell{1},n1,n2,n3,idxnz,1);
        X0{2} = Wavefun(Qcell{1},n1,n2,n3,idxnz,2);
    end
end

end
