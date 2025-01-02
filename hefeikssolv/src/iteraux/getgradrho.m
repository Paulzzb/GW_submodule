function [grho,grho2] = getgradrho(mol,rho,icomplex)
% GETGRADRHO gradient of rho.
%    [grho,grho2] = GETGRADRHO(mol,rho) returns the gradient and gradient
%    absolute square of the grho.
%
%   See also getcharge.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if (nargin == 2)
  icomplex = 0;
end

n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;
ecut2 = mol.ecut2;

grid2 = Ggrid(mol,ecut2);
idxnz2 = grid2.idxnz;
gkx2 = grid2.gkx;
gky2 = grid2.gky;
gkz2 = grid2.gkz;

% Calculate Hartree potential
rhog   = fft3(rho);
rhog   = rhog(idxnz2);

grhog = zeros(n1*n2*n3,3);
grhog(idxnz2,1) = 1i*rhog.*gkx2;
grhog(idxnz2,2) = 1i*rhog.*gky2;
grhog(idxnz2,3) = 1i*rhog.*gkz2;
if icomplex == 0
  grho = real(ifft3(reshape((grhog),n1,n2,n3,3)));
else
  grho = ifft3(reshape((grhog),n1,n2,n3,3));
end

grho2 = sum(grho.^2,4);

end
