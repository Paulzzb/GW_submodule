function vhart = getVhart(mol,rho)
% GETVHART Calculates the Hartree potential.
%    vhart = getVhart(mol,rho) calculates the Hartree potential
%    induced by electron repulsion interation.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% extracting the discretiation dimension in real space
n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;
% energy cut-off from the molecule
ecut2 = mol.ecut2;
% add temperarily
grid2  = Ggrid(mol,ecut2);
% indeces of the non-zero elements in Fourier space
idxnz2 = grid2.idxnz;
% Fourier grid
gkk2   = grid2.gkk;

if mol.nspin == 1
    arho = rho;
elseif mol.nspin == 2
    arho = rho{1}+rho{2};
elseif mol.nspin == 4
    arho = rho{1};
end

% Calculate Hartree potential
rhog   = fft3(abs(arho));
rhog   = rhog(idxnz2);

w      = zeros(n1,n2,n3);
inz    = abs(gkk2) ~= 0;
w(idxnz2(inz)) = e2Def()*4*pi*rhog(inz)./gkk2(inz);
vhart  = real(ifft3(w));

end
