function vtot = getVtot(mol,vion,vext,vhart,vxc)
%
% Usage: vtot = getvtot(mol,vion,vext,vhart,vxc);
%
% Purpose:
%   Compute the total potential energy and keep it band limited
%
% Input: 
%     mol   a Molecular object
%    vion   a 3D array that contains the local ionic potential
%    vext   a 3D array that contains external potential. In most
%           case, this should be an empty array
%    vhart  a 3D array that contains the Hartree (Coulomb) potential.
%    vxc    a 3D array or cell array that contains the exchange-correlation potential
%
% Output:
%    vtot   a 3D array or cell array that contains the total (local) potential
%
nspin = mol.nspin;
[n1,n2,n3]=size(vion);
ecut2   = mol.ecut2;
grid2  = Ggrid(mol,ecut2);
idxnz2  = grid2.idxnz;
if nspin == 1||(mol.nspin == 4 && ~mol.domag)
    vtot = vion + vhart + vxc;
    if (~isempty(vext))
      vtot = vtot + vext;
    end
elseif nspin == 2
    vtot = cell(2,1);
    for is = 1:2
        vtot{is} = vion + vhart + vxc{is};
        if (~isempty(vext))
            vtot{is} = vtot{is} + vext;
        end
    end
elseif nspin == 4 && mol.domag
    vtot = cell(4,1);
    vtot{1} = vion + vhart + vxc{1};
    for i = 2:4
        vtot{i} = vxc{i};
    end
    if ~isempty(vext)
        for i = 1:4
            vtot{i} = vtot{i} + vext;
        end
    end
end
% keep vtot band limited (sinc interpolation)
gm2 = zeros(n1,n2,n3);
gm2(idxnz2) = 1;
iz = gm2==0;
if nspin == 1
    vfft = fftn(vtot);
    vfft(iz) = 0;
    vtot = ifftn(vfft);
elseif nspin == 2||nspin == 4
%     for is = 1:nspin
%         vfft = fftn(vtot{is});
%         vfft(iz) = 0;
%         vtot{is} = ifftn(vfft);   
%     end    
end
