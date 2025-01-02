function [diagE_ISDF] = tddft_couplings_isdf(ksinfo, mol, nvbands, ncbands, vcrank_ratio)
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

timeforFXC = toc

% size of psivr: nr * nvbands; size of psicr: nr * ncbands
psivr = zeros(nr,nvbands);
psicr = zeros(nr,ncbands);
psivr = F'*Z(:,nv-nvbands+1:nv);
psicr = F'*Z(:,nv+1:nv+ncbands);

% size of vcrzeta_mu: nr * n_mu
% size of psivc_mu: n_mu * nvc (nc *nv)
% size of vczeta_mu: ng * n_mu
vcrank_mu = ceil(nvbands*ncbands*vcrank_ratio);
n_mu = vcrank_mu;

vcrzeta_mu = zeros(nr, vcrank_mu);
vczeta_mu = zeros(ng, vcrank_mu);
psivc_mu = zeros(ng, vcrank_mu);

vcrzeta_mu = zeros(nr, vcrank_mu);

%timeISDF = clock;
%[vcrzeta_mu, vcind_mu] = isdf(psicr, conj(psivr), vcrank_mu);
%timeforISDF = etime(clock,timeISDF)

timeISDF = clock;
[vcrzeta_mu, vcind_mu] = isdf_fast(psicr, conj(psivr), vcrank_mu);
timeforISDF = etime(clock,timeISDF)

tic

vczeta_mu = F*vcrzeta_mu;

timeforFFT = toc

tic

psivc_mu = prod_states(psicr(vcind_mu, :), conj(psivr(vcind_mu, :)));

timeforCmu = toc

n_mu = vcrank_mu;


if(1)

% Coulomb couplings 
fac = 2.0*vol;
nvc = nvbands*ncbands;
A1 = zeros(vcrank_mu,vcrank_mu);
Dcoulg = spdiags(ksinfo.coulG,0,ng,ng);

tic

A1 = vczeta_mu'*Dcoulg*vczeta_mu*fac;

timeforCoulomb = toc

% Exchange-correlation couplings 
[m1,m2,m3] = size(fxc);
m123 = m1*m2*m3;
fxc2 = reshape(fxc,m123,1);
fxc2r = spdiags(fxc2,0,nr,nr);
fac = 2.0*vol.^3/m123;
A2 = zeros(vcrank_mu,vcrank_mu);

tic

A2 = vcrzeta_mu'*fxc2r*vcrzeta_mu*fac;

timeforXC = toc

end %if(1)

if(0)

% Coulomb couplings 
fac = 2.0*vol;
nvc = nvbands*ncbands;
A1 = zeros(vcrank_mu,vcrank_mu);

tic

for j1 = 1 : vcrank_mu
  for j2 = 1 : vcrank_mu
    A1(j1, j2) = sum(conj(vczeta_mu(:,j1)).*(Dcoul.*vczeta_mu(:,j2))) * fac;
  end
end

timeforCoulomb = toc

% Exchange-correlation couplings 
[m1,m2,m3] = size(fxc);
m123 = m1*m2*m3;
fxc2 = reshape(fxc,m123,1);
fac = 2.0*vol.^3/m123;
A2 = zeros(vcrank_mu,vcrank_mu);

tic

for j1 = 1 : vcrank_mu
  for j2 = 1 : vcrank_mu
    A2(j1, j2) = sum(conj(vcrzeta_mu(:,j1)).*(fxc2.*vcrzeta_mu(:,j2))) * fac;
  end
end

timeforXC = toc

end %if(0)

A = zeros(nvc,nvc);

tic

A = psivc_mu'*(A1+A2)*psivc_mu;

timeforC = toc

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

diagE = eig(A);
diagE = 2.0*diagE;

timeforDiagHxc = toc

%diagE2 = reshape(diagE,[ncbands,nvbands]);

tic

diagE_ISDF = sort(diagE);

timeforsortdiagE = toc
