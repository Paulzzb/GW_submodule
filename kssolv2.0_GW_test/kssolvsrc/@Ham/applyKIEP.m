function hpsi = applyKIEP(H,X)
%
% Perform the multiplication of the Kinetic and Ionic
% potential part of the KS Hamiltonian with
% wavefunction(s) X.
%
% usage: KX = applyKIEP(H,X);
%
n1 = X.n1;
n2 = X.n2;
n3 = X.n3;

gkin = H.gkin;
idxnz  = H.idxnz;

psir = ifft3(X);
psir = bsxfun(@times,H.vion(:)+H.vext(:),psir);
psir = fft3(psir);
if X.iscompact
    psir = compress(psir,idxnz);
else
    gm = ones(n1*n2*n3,1);
    gm(idxnz) = 0;
    psir(gm,:) = 0;
end

hpsi = psir + bsxfun(@times,gkin,X);
vnlmat  = H.vnlmat;
vnlsign = H.vnlsign;
hpsi = hpsi + vnlmat*bsxfun(@times,vnlsign,(vnlmat'*X));

end