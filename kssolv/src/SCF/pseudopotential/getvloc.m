function [vion,rho] = getvloc(mol)
% GETVLOC returns the local pseudo potential of the molecule.
%    [vion,rho] = GETVLOC(mol) calculate the local pseudo potential of the
%    molecule based on the local pseudo potential of each atom in the
%    molecule in G-space. And return the corresponding atomic charge
%    density as well.
%
%    See also vnl2g, vloc2g, getvnl. 

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

n1       = mol.n1;
n2       = mol.n2;
n3       = mol.n3;
vol      = mol.vol;
atoms    = mol.atoms;
nspin   = mol.nspin;

vion  = zeros(n1,n2,n3);
if nspin == 1
    rho   = zeros(n1,n2,n3);
else
    rho = cell(nspin,1);
    for i = 1:nspin
        rho{i} = zeros(n1,n2,n3);
    end
end

if mol.domag
    mag = mol.initmag;
end

if isempty(atoms)
    return;
end

xyzlist  = mol.xyzlist;
ntypes = length(mol.natoms);
alist = mol.alist;

vloc3dg = zeros(n1,n2,n3);
if nspin == 1
    rho3dg  = zeros(n1,n2,n3);
else
    rho3dg = cell(nspin,1);
    for i = 1:nspin
        rho3dg{i} = zeros(n1,n2,n3);
    end
end
ppvar = mol.ppvar;
vlocg = ppvar.vlocg;
rhog  = ppvar.rhog;

ecut2  = mol.ecut2;
grid2 = Ggrid(mol,ecut2);
inz    = grid2.idxnz;
gkx2   = grid2.gkx;
gky2   = grid2.gky;
gkz2   = grid2.gkz;

% vectorize the original for i = 1:nnz loop
for itype = 1:ntypes
    index  = alist == itype;
    xyzmat = xyzlist(index,:)';
    phmat  = [gkx2 gky2 gkz2]*xyzmat;    
    % ccvec is the structure factor used to account for
    % the translation of the atomic position    
    ccvec  = exp(-1i.*phmat);
    vloc3dg(inz)  = vloc3dg(inz) + vlocg(:,itype).*sum(ccvec,2);
    if nspin == 1
        rho3dg(inz)  = rho3dg(inz) + rhog(:,itype).*sum(ccvec,2);
    elseif nspin == 2
        rho3dg{1}(inz) = rho3dg{1}(inz) + rhog(:,itype).*sum(ccvec,2);       
        % rho{1} is total density read from pseudopotential files
        % rho{2} is determined by initial magnetization.
        amag = repmat(mag.amag(index),length(inz),1);
        rho3dg{2}(inz) = rho3dg{2}(inz) + rhog(:,itype).*sum(ccvec.*amag,2);
    elseif nspin == 4
        rho3dg{1}(inz) = rho3dg{1}(inz) + rhog(:,itype).*sum(ccvec,2);
        % rho{1} is total density read from pseudopotential files
        % rho{2~4} is determined by initial magnetization.
        if mol.domag
            amag = repmat(mag.amag(index),length(inz),1);
            angular = cell(3,1);
            angular{1} = sin(mag.theta(index)).*cos(mag.phi(index));
            angular{2} = sin(mag.theta(index)).*sin(mag.phi(index));
            angular{3} = cos(mag.theta(index));
            for is = 2:nspin
                rho3dg{is}(inz) = rho3dg{is}(inz) + rhog(:,itype).*sum(ccvec.*amag...
                    .*repmat(angular{is-1},length(inz),1),2);
            end
        end
    end
end

% get rho in real space 
vion  = real(ifft3(vloc3dg)*n1*n2*n3)*4/vol*pi;
if nspin == 1
    rho   = real(ifft3(rho3dg)*n1*n2*n3)/vol;
elseif nspin == 2
    for i = 1:nspin
        rho_real  = real(ifft3(rho3dg{i})*n1*n2*n3)/vol;
        rho{i}     = rho_real;
    end
    rho = transform_rho(mol,rho,'sum2updw');
elseif nspin == 4
    for i = 1:nspin
        rho_real  = real(ifft3(rho3dg{i})*n1*n2*n3)/vol;
        rho{i}     = rho_real;
    end
end
