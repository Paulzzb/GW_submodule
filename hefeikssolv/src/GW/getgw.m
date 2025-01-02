function [Esx_x, Ech, Ex] = getgw(ksinfo, eta, nvbands, ncbands, chi0,  QE_WAV)
%
%
% ksinfo is a structure that contains ground state calculation results
% and parameters of the molecules
%   nv      --- number of valence states
%   Z       --- contains the eigenvecgtors from the KSDFT calculation
%   ev      --- contains the corresponding eigenvalues
%   vol     --- volume of the unit cell
%   ntot    --- total number of grid points on which the wavefunction is sampled
% 
% eta     --- Lorentzian broadening factor that turns a Dirac delta into a 
%             a smoother peak
% nvbands --- the number of valence bands (KS orbitals) from the Fermi 
%             level included in the kernel calculation. 
% ncbands --- the number of empty bands (KS orbitals) included in 
%             the kernel calculation
%
%
%
format long;

nv   = ksinfo.nv;
if nargin>5
Z    = QE_WAV;
else
Z    = ksinfo.Z;
end
ev   = ksinfo.ev;
F    = ksinfo.F;
vol  = ksinfo.vol;
ntot = ksinfo.ntot;
%
% 1:nv = occupied, nv+1:n = unoccupied
im = sqrt(-1);
[ng,nr]=size(F);

nbands = nvbands + ncbands;

Z'*Z;
% Normalize the wavefunction in the Fourier space
% 1/Vol sum_G |Psi(G)|^2 = 1
for iv = 1 : nbands
  Z(:,iv) = Z(:,iv) * sqrt (vol) / (norm(Z(:,iv)));
end

ng = size(chi0,1);
I = eye(ng);
vol  = ksinfo.vol;

Dcoul0 = spdiags(ksinfo.coulG(:,4),0,ng,ng);

epsilon = geteps(chi0, ksinfo.coulG);

inveps = inv(epsilon);
I = eye(ng);
inveps_1 = inveps-I;    % epsilon^-1 - I
ksinfo.coulG(1,4) = ksinfo.coulG0;

Dcoul   = spdiags(ksinfo.coulG(:,4),0,ng,ng);
W1 = inveps_1*Dcoul;   %calculate SX-X

Esx_x = zeros(nbands, nbands);
Ech = zeros(nbands, 1);
Ex = zeros(nbands, nbands);

for n1 = 1:nbands
  psin1r = conj(F'*Z(:,n1));
  for n2 = 1:nbands
    psin2r = conj(F'*Z(:,n2));
    for iocc = 1:nv
      psioccr = F'*Z(:,iocc);
      psioccn1 = F*(psioccr.*psin1r);
      psioccn2 = F*(psioccr.*psin2r);
      %Esx(n1, n2) = Esx(n1, n2) + psioccn1'*inveps*Dcoul*psioccn2;
      Esx_x(n1, n2) = Esx_x(n1, n2) + psioccn1'*(W1*psioccn2);
      Ex(n1, n2) = Ex(n1, n2) + psioccn1'*(Dcoul*psioccn2);
    end  
  end
end

%%%%%calculate CH  %%%%%%%%%%

%calculate aqsch
for iv1 = 1:nbands
  psiiv1 = F'*Z(:,iv1);
  psiiv2 = F'*Z(:,iv1);
  aqsch(:,iv1) = F*(conj(psiiv1).*psiiv2);
end

for m=1:nbands
  for i=1:ng
  schx_tmp = 0;
    for j=1:ng
      if abs(inveps_1(i,j)) > 1E-6
        G=ksinfo.coulG(i,1:3)-ksinfo.coulG(j,1:3);
        for k=1:ng
          if G==ksinfo.coulG(k,1:3)
                schx = aqsch(k,m) * inveps_1(j,i);
                schx_tmp = schx_tmp + schx;
          end 
        end
      end
    end
    Ech(m) = Ech(m) + schx_tmp * ksinfo.coulG(i,4);
  end
end

fac = 1.00/vol;

Esx_x = 0.0 - fac * Esx_x;
Ech = 0.5 * fac * Ech;
Ex = 0.0 - fac * Ex;
