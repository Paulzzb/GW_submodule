function HX = mtimes(H,X)
% MTIMES  Overload multiplication operator for Ham class
%    HX = H*X returns a wavefun corresponding to H*X.
%
%    See also Ham, Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

ncol = ncols(X);
npw = numel(X.idxnz);
npol = size(X.psi,1)/npw;
nspin = H.nspin;
lspinorb = H.lspinorb;

% Apply Laplacian
KinX = repmat(H.gkin,npol,ncol).*X;

% Apply total local potential
if npol == 1
    Xr3d  = ifft3(X);
    if nspin == 1
        vXr3d = repmat(H.vtot(:),1,ncol).*Xr3d;
    elseif nspin == 2
        vtot = H.vtot{X.ispin};
        vXr3d = repmat(vtot(:),1,ncol).*Xr3d;
    end
    VX    = fft3(vXr3d);
    VtotX = VX(X.idxnz,:);
elseif npol == 2
    Xup = Wavefun(X.psi(1:npw,:),X.n1,X.n2,X.n3,X.idxnz);
    Xdw = Wavefun(X.psi(npw+1:end,:),X.n1,X.n2,X.n3,X.idxnz);
    Xrup = ifft3(Xup);
    Xrdw = ifft3(Xdw);
    if iscell(H.vtot)
        sup = repmat(H.vtot{1}(:)+H.vtot{4}(:),1,ncol).*Xrup + ...
            repmat(H.vtot{2}(:)-1i*H.vtot{3}(:),1,ncol).*Xrdw;
        sdw = repmat(H.vtot{1}(:)-H.vtot{4}(:),1,ncol).*Xrdw + ...
            repmat(H.vtot{2}(:)+1i*H.vtot{3}(:),1,ncol).*Xrup;
    else
        sup = repmat(H.vtot(:),1,ncol).*Xrup;
        sdw = repmat(H.vtot(:),1,ncol).*Xrdw;
    end
    VXup = fft3(sup);
    VXdw = fft3(sdw);
    VtotX = zeros(npw*npol,ncol);
    VtotX(1:npw,:) = VXup.psi(X.idxnz,:);
    VtotX(npw+1:end,:) = VXdw.psi(X.idxnz,:);
end
% Apply nonlocal pseudopotential
if npol == 1
    VnlX = H.vnlmat*(repmat(H.vnlsign,1,ncol).*(H.vnlmat'*X));
elseif npol == 2
    becp_up = H.vnlmat'*Xup;
    becp_dw = H.vnlmat'*Xdw;

    if lspinorb
        psup = H.vnlsign{1}*becp_up + H.vnlsign{2}*becp_dw;
        psdw = H.vnlsign{3}*becp_up + H.vnlsign{4}*becp_dw;
    else
        % H.vnlsign{2}{3} = 0 without soc
        psup = H.vnlsign{1}*becp_up;
        psdw = H.vnlsign{4}*becp_dw;
    end
    VnlX = zeros(npw*npol,ncol);
    VnlX(1:npw,:) = H.vnlmat*psup;
    VnlX(npw+1:end,:) = H.vnlmat*psdw;
end

HX = KinX + VtotX + VnlX;
% Apply Fock exchange operator for hybrid functional
if H.ishybrid
    if ~isempty(X.ik) && X.ik~=H.ik
        error('k points of H and X are different! %d %d',X.ik,H.ik);
    end
    X.ik = H.ik;
    X.wks= H.wks;
    Vexx = H.vexx;
    HX = HX + Vexx(X);
end

end
