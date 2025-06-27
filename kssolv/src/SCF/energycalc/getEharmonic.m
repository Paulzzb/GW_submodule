function Eharmonic = getEharmonic(rho_in, vHarmonic)
% GETEHARMONIC Computes a Harmonic potential energy.
%    EHarmonic = getEHarmonic(rho_in, vHarmonic) get the Harmonic potential energy
%    by EHarmonic = int dr V(r) * rho_in(r) = (1/2) w**2 * int dr { (r-r0)**2 rho_in(r) }
%    where r0 is the center of the cell.
%
%    rho_i     --- input charge density
%    vHarmonic --- Harmonic potential (from getvHarmonic) 

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

[n1,n2,n3] = size(rho_in);
[m1,m2,m3] = size(vHarmonic);
if (m1 ~= n1 | m2 ~= n2 | m3 ~= n3)
   error('getEHarmonic: Dimension of rho_in and vHarmonic do not match.');
end

n123 = n1*n2*n3;

rho1d = reshape(rho_in,n123,1);
v1d = reshape(vHarmonic,n123,1);
Eharmonic = rho1d'*v1d;   % Check the normalization.

end
