function [tab,D] = vnl2g(pp,ecut,noncolin)
% VNL2G converts non-local pseudo potential at irregular radial grid in
%    real-space to uniform radial grid in G-space.
%    [tab,D] = VNL2G(pp,ecut) calculate the Fourier transform of the
%    non-local pseudo potential to G-space at uniform radial grid in the
%    G-space. The input of this function is the pseudo potential and energy
%    cut. And the function returns the Vnl at G-space in tab and the
%    corresponding middle matrix D. The detailed implementation of the
%    Fourier transform are detailed in the documents.
%
%    See also VLOC2G, GETWQ.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

dq = 0.01; % gap for interpolation
nqxq = round(sqrt(ecut*2*meDef())/dq+4);

nb = pp.nonloc.nbeta;
beta = pp.nonloc.beta;
hbeta = min(max(pp.nonloc.cutoff_radius_index),length(pp.r));
r = pp.r(1:hbeta);
rab = pp.rab(1:hbeta);
sD = pp.nonloc.dij;
lll = pp.nonloc.lll;

lenl = length(lll);
totall = sum(2*lll+1);

tab = zeros(nqxq,nb);

if nb ==0
    D = [];
    return;
end

lmaxx = 4;
lqmax = 2*lmaxx+1;
rot_ylm = zeros(lqmax,lqmax);

for it = 1:nb
    besr = bessel(pp.nonloc.lll(it),r*dq*((1:nqxq)-1));
    aux = repmat(beta(1:hbeta,it),1,nqxq).*besr.*repmat(r,1,nqxq);
    tab(:,it) = simpson(hbeta,aux,rab);
end 

% In the spin-orbit case we need the unitary matrix u which rotates the 
% real spherical harmonics and yields the complex one.
if pp.info.has_so
    l = lmaxx;
    rot_ylm(l+1,1) = complex(1,0);
    for n1 = 2:2:2*l+1
        m = n1/2;
        n = l+1-m;
        rot_ylm(n,n1) = complex((-1)^m/sqrt(2),0);
        rot_ylm(n,n1+1) = complex(0,-(-1)^m/sqrt(2));
        n = l+1+m;
        rot_ylm(n,n1) = complex(1/sqrt(2),0);
        rot_ylm(n,n1+1) = complex(0,1/sqrt(2));
    end
end

% first calculate the fcoef coefficients
if pp.info.has_so
    fcoef = zeros(totall,totall,2,2);
    iidx = 0;
    for iit = 1:lenl
        ill = pp.nonloc.lll(iit);
        ijj = pp.spinorb.jjj(iit);
        for iitm = -ill:ill
            iidx = iidx+1;
            jidx = 0;
            for jit = 1:lenl
                jll = pp.nonloc.lll(jit);
                jjj = pp.spinorb.jjj(jit);
                for jitm = -jll:jll
                    jidx = jidx+1;
                    if ill==jll && ijj==jjj
                        for is1 = 1:2
                            for is2 = 1:2
                                coeff = 0;
                                for m = -ill-1:ill
                                    m0 = sph_ind(ill,ijj,m,is1)+lmaxx+1;
                                    m1 = sph_ind(jll,jjj,m,is2)+lmaxx+1;
                                    coeff = coeff+rot_ylm(m0,iitm+ill+1)*spinor(ill,ijj,m,is1)...
                                        *conj(rot_ylm(m1,jitm+jll+1))*spinor(jll,jjj,m,is2);
                                end
                                fcoef(iidx,jidx,is1,is2) = coeff;
                            end
                        end
                    end
                end
            end
        end
    end
    % then calculate the bare coefficients
    D = zeros(totall,totall,4);
    iidx = 0;
    for iit = 1:lenl
        ill = pp.nonloc.lll(iit);
        for iitm = -ill:ill
            iidx = iidx+1;
            jidx = 0;
            for jit = 1:lenl
                jll = pp.nonloc.lll(jit);
                for jitm = -jll:jll
                    jidx = jidx+1;
                    ijs = 0;
                    for is1 = 1:2
                        for is2 = 1:2
                            ijs = ijs+1;
                            D(iidx,jidx,ijs) = sD(iit,jit)*fcoef(iidx,jidx,is1,is2);
                            if iit == jit
                                fcoef(iidx,jidx,is1,is2) = 0;
                            end
                        end
                    end
                end
            end
        end
    end
else
    % TODO: The following loop could be reduced. Since nb and lll are all
    % small, it might be unnecessary to recude it.
    if ~noncolin
        D = zeros(totall);
    else
        D = zeros(totall,totall,2,2);
    end
    iidx = 0;
    for iit = 1:lenl
        ill = pp.nonloc.lll(iit);
        for iitm = -ill:ill
        iidx = iidx+1;
        jidx = 0;
            for jit = 1:lenl
                jll = pp.nonloc.lll(jit);
                for jitm = -jll:jll
                    jidx = jidx+1;
                    if ill==jll && iitm==jitm
                        if ~noncolin
                            D(iidx,jidx) = sD(iit,jit);
                        else 
                            % as for noncolinear-spin without soc, only first and last D are nonzero
                            D(iidx,jidx,1) = sD(iit,jit);
                            D(iidx,jidx,4) = sD(iit,jit);
                        end
                    end
                end
            end
        end
    end   
end
end
