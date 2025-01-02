function [diagE] = tddft_couplings_gpu(ksinfo, mol, nvbands, ncbands)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization 
nv      = ksinfo.nv;
Z       = ksinfo.Z;
Z       = gpuArray(Z);
ev      = ksinfo.ev;
F       = ksinfo.F;
vol     = ksinfo.vol;
ntot    = ksinfo.ntot;
rho     = ksinfo.rho;
[ng,nr] = size(F);

% Normalize the wavefunction in the Fourier space
% sum_G |Psi(G)|^2 = 1

tic

for iv = 1 : nv+ncbands
  Z(:,iv) = Z(:,iv) / (norm(Z(:,iv)));
end

timeforZ=toc

% Coulomb kernel in reciprocal space
Dcoul = ksinfo.coulG;

% Exchange-correlation kernel in real space

tic

fxc = getfxc(mol,rho);

fxc = gpuArray(fxc);

timeforFXC = toc

if(1)
   
psivcg = gpuArray.zeros(ng, nvbands*ncbands);
psivcr = gpuArray.zeros(nr, nvbands*ncbands);

tic  
for iv = 1:nvbands
  psivr = conj(F'*Z(:,nv-nvbands+iv));
  for ic = 1:ncbands
    ivc = (iv-1)*ncbands+ic;
    psivcr(:,ivc) = conj((F'*Z(:,nv+ic)).*psivr);
    psivcg(:,ivc) = F*psivcr(:,ivc);
  end
end

timeforMij = toc

fac = 2.0*vol;
nvc = nvbands*ncbands;
A1 = gpuArray.zeros(nvc,nvc);
Dcoulg = spdiags(ksinfo.coulG,0,ng,ng);
Dcoulg = gpuArray(Dcoulg);

tic
  
A1 = psivcg'*(Dcoulg*psivcg)*fac;

timeforCoulomb = toc

[m1,m2,m3] = size(fxc);
m123 = m1*m2*m3;

tic

fxc2 = reshape(fxc,m123,1);

timeforFXC2 = toc

fac = 2.0*vol.^3/m123;
fxc2r = spdiags(fxc2,0,nr,nr);
A2 = gpuArray.zeros(nvc,nvc);

tic

A2 = psivcr'*(fxc2r*psivcr)*fac;


timeforXC = toc

A = gpuArray.zeros(nvc,nvc);
A = A1 + A2;

end %if(1)

if(0)

% Coulomb couplings 
fac = 2.0*vol;
nvc = nvbands*ncbands;
A = zeros(nvc,nvc);

tic

for n1 = 1 : nvbands
  psin1r = conj(F'*Z(:,nv-nvbands+n1)); % n123/vol
  for n2 = 1 : ncbands
    psiocc12 = F*(F'*Z(:,nv+n2).*psin1r); % (vol/n123) * (n123/vol)^2
    for n3 = 1 : nvbands
      psin3r = conj(F'*Z(:,nv-nvbands+n3));
      for n4 = 1 : ncbands
        psiocc34 = F*(F'*Z(:,nv+n4).*psin3r);
        n12 = (n1-1)*ncbands+n2;
        n34 = (n3-1)*ncbands+n4;
        A(n12, n34) = A(n12, n34) + sum(conj(psiocc34).*(Dcoul.*psiocc12)) * fac;
      end
    end  
  end
end

timeforCoulomb = toc

% Exchange-correlation couplings 
[m1,m2,m3] = size(fxc);
m123 = m1*m2*m3;
fxc2 = reshape(fxc,m123,1);
fac = 2.0*vol.^3/m123;

tic

for n1 = 1 : nvbands
  psin1r = conj(F'*Z(:,nv-nvbands+n1)); % n123/vol
  for n2 = 1 : ncbands
    psiocc12 = F'*Z(:,nv+n2).*psin1r;   % (n123/vol)^2
    for n3 = 1 : nvbands
      psin3r = conj(F'*Z(:,nv-nvbands+n3));
      for n4 = 1 : ncbands
        psiocc34 = F'*Z(:,nv+n4).*psin3r;
        n12 = (n1-1)*ncbands+n2;
        n34 = (n3-1)*ncbands+n4;
        A(n12, n34) = A(n12, n34) + sum(conj(psiocc34).*(fxc2.*psiocc12)) * fac;
      end
    end  
  end
end

timeforXC = toc

end %if(0)

tic

% Diagonal part
for iv = 1 : nvbands
  for ic = 1 : ncbands
    ivc = (iv-1)*ncbands+ic;
    A(ivc,ivc) = A(ivc,ivc) + ev(nv+ic) - ev(nv-nvbands+iv);
  end
end

timeforHxc = toc

% Diagonalization

tic
A     = gpuArray(A);
diagE = eig(A);
diagE = 2.0*diagE;
A     = gather(A);

timeforDiagHxc = toc

%tic
%ncol = min(nvbands*4,nvbands*ncbands);
%diagE = eigs(A,ncol,'SR'); 
%diagE = 2.0*diagE;
%timeforDiagHxc2 = toc

%diagE2 = reshape(diagE,[ncbands,nvbands]);

tic

diagE = sort(diagE);

timeforsortdiagE = toc
