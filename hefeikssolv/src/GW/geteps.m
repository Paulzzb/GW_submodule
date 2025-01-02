function epsilon = geteps(chi0, coulG)
%
%  epsilon = geteps(chi0, coulG);
% 
%  compute the dieletric matrix epsilon from irreducible polarizability chi0 
%  and the Coulomb interaction in reciprocal space
%
ng = size(chi0,1);
I = eye(ng);
%coulG(1,4) = 314.15926;
Dcoul = spdiags(coulG(:,4),0,ng,ng);
epsilon = I - Dcoul*chi0;
