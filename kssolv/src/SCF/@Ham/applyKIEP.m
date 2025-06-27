function hpsi = applyKIEP(H,X)
% APPLYKIEP Perform the multiplication of the Kinetic and Ionic
% potential part of the KS Hamiltonian with wavefunction(s) X.
%    hpsi = applyKIEP(H,X) returns the Kinetic and Ionic potential 
%    applied wavefunction
%
% 	 See also Ham, Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

n1 = X.n1;
n2 = X.n2;
n3 = X.n3;
ncol = ncols(X);

gkin = H.gkin;
idxnz  = H.idxnz;

psir = ifft3(X);
vpsir = repmat(H.vion(:)+H.vext(:),1,ncol).*psir;
vpsi = fft3(vpsir);
if X.iscompact
    vpsi = compress(vpsi,idxnz);
else
    gm = ones(n1*n2*n3,1);
    gm(idxnz) = 0;
    vpsi(gm,:) = 0;
end

hpsi = vpsi + repmat(gkin,1,ncol).*X;
vnlmat  = H.vnlmat;
vnlsign = H.vnlsign;
hpsi = hpsi + vnlmat*(repmat(vnlsign,1,ncol).*(vnlmat'*X));

end
