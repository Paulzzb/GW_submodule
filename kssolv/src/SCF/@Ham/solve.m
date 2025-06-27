function X = solve(H,B,shift)
% SOLVE Function to solve (H - shift*I)X = B
% where H is a Hamiltonian
% 	X = solve(H,B,shift) returns the solution of (H - shift*I)X = B
%
%   See also Ham, Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.	

if nvargin == 2
    shift = zeros(1,ncol(B));
end

n1 = H.n1;
n2 = H.n2;
n3 = H.n3;
idxnz = H.idxnz;

xArray = zeros(length(idxnz), ncol(B));

for ii=1:ncol(B)
    % solving the problems using GMRES
    xArray(:,ii) = gmres( @(p) shift(ii)*p-mult(H,p), B.psi, [], 1e-6, 300 );
    
end
% Building the resultign Wavefield
X = Wavefun(xArray,n1,n2,n3,idxnz);

end
