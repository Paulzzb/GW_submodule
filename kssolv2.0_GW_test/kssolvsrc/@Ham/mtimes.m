function HX = mtimes(H,X)
% MTIMES  Overload multiplication operator for Ham class
%    HX = H*X returns a wavefun corresponding to H*X.
%
%    See also Ham, Wavefun.

% Apply Laplacian
KinX = bsxfun(@times,H.gkin,X);

% Apply total local potentials
VtotX = ifft3(X);
VtotX = bsxfun(@times,H.vtot(:),VtotX);
VtotX = fft3(VtotX);
VtotX = compress(VtotX,X.idxnz);

% Apply nonlocal pseudopotential
VnlX = H.vnlmat*bsxfun(@times,H.vnlsign,(H.vnlmat'*X));

HX = KinX + VtotX + VnlX;

if H.ishybrid
  Vexx = H.vexx;
  HX = HX + Vexx(X);
end
end
