function vtot = getVtot(vion,vext,vhart,vxc)
%
% Usage: vtot = getvtot(vion,vext,vhart,vxc);
%
% Purpose:
%   Compute the total potential energy and keep it band limited
%
% Input: 
%    vion   a 3D array that contains the local ionic potential
%    vext   a 3D array that contains external potential. In most
%           case, this should be an empty array
%    vhart  a 3D array that contains the Hartree (Coulomb) potential.
%    vxc    a 3D array that contains the exchange-correlation potential
%
% Output:
%    vtot   a 3D array that contains the total (local) potential
%

vtot = vion + vhart + vxc;
if (~isempty(vext))
  vtot = vtot + vext;
end

end
