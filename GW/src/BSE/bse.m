function [d,Vs] = bse(GWinfo,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%construct Hamiltonian of BSE, solving the BSE eigenvalue and eigenvector,
%calculate the absorption without and with electron-hole interaction.
%
%GWinfo is a structure that contains ground state calculation results
%and parameters of the molecules.
%   nv           --- number of valence states
%   Z            --- contains the eigenvecgtors from the KSDFT calculation or Quantum-Espresso calculation
%   ev           --- contains the corresponding eigenvalues
%   vol          --- volume of the unit cell
%   ntot         --- total number of grid points on which the wavefunction is sampled
%   F            --- Fourier tranform
%   Vxc          --- exchange and correlation of DFT
%   eqp          --- quasiparticle energy
%   bdot         --- reciprocal lattice
%   G            --- the grid in reciprocal lattice
%   plasma_omega --- Plasma Frequency 4.0*sqrt(ne*pi/vol)
%   eta_abs      --- Lorentzian broadening factor that turns a Dirac delta into a 
%                    a smoother peak
%   eta          --- Lorentzian broadening factor that turns a Dirac delta into a 
%                    a smoother peak
%   nvbands      --- the number of valence bands (KS orbitals) from the Fermi 
%                    level included in the kernel calculation. 
%   ncbands      --- the number of empty bands (KS orbitals) included in 
%                    the kernel calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long g;

%maxNumCompThreads(1);

% Initialization
nv    = GWinfo.nv;
Z     = GWinfo.Z;
ev    = GWinfo.ev;
F     = GWinfo.F;
vol   = GWinfo.vol;
Vxc   = GWinfo.Vxc;
eqp   = GWinfo.eqp;
bdot  = GWinfo.bdot;
G     = GWinfo.G(:,1:3);
plasma_omega = GWinfo.plasma_omega;
nvbands = GWinfo.nvbands;
ncbands = GWinfo.ncbands;
nbands = nvbands + ncbands;
NBANDS = GWinfo.NBANDS;
ry2ev = 13.60569253;

%parameters for head/wing/body/exchange of bse Hamiltonian
fac_d = - 8*pi/vol;
w_eff = GWinfo.coulG0/(8*pi);
hbse_mat_head = fac_d*w_eff;
hbse_mat_wing = fac_d;
hbse_mat_body = fac_d;
hbse_mat_exchange = -fac_d*2;

%F = KSFFT(mol);
im = sqrt(-1);
[ng,nr] = size(F);

fprintf('nr = %d, ng = %d, nbands = %d\n', nr, ng, nbands);

% Normalize the wavefunction in the Fourier space
% 1/Vol sum_G |Psi(G)|^2 = 1

%startGW = tic;

startZ = tic;

if NBANDS == 0
  for iv = nv-nvbands+1:nbands
    Z(:,iv) = Z(:,iv) * sqrt(vol) / (norm(Z(:,iv)));
  end
else
  for iv = nv-nvbands+1:NBANDS
    Z(:,iv) = Z(:,iv) * sqrt(vol) / (norm(Z(:,iv)));
  end
end

timeforZ = toc(startZ)

if(1)

startC2R = tic;

if NBANDS == 0
  psir = zeros(ng,nbands);
  psir = F'*Z(:,nv-nvbands+1:nbands);
else 
  psir = zeros(ng,NBANDS);
  psir = F'*Z(:,nv-nvbands+1:NBANDS);
end

timeforC2R = toc(startC2R)

chi0 = zeros(ng);

startchi0 = tic;

psivcr = prod_states_gw((psir(:,nv-nvbands+1:nv)),psir(:,nv+1:nbands),ev,nv);
Mg = F*psivcr;
chi0 =  4.0 * Mg*Mg'/vol;

timeforchi0 = toc(startchi0)

startW = tic;

I = eye(ng);
%chi0(logical(eye(ng)))=diag(chi0).*GWinfo.coulG(:,4);
%chi0_coul=chi0;
%inveps = inv(I + chi0_coul)
%chi0.*GWinfo.coulG(:,4)
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
%Dcoul*chi0
inveps = inv(I + Dcoul*chi0);
%inveps-I
%inveps = epsilon\I;
%diag(inveps)

W = inveps * Dcoul/8.0/pi;

timeforW = toc(startW)

startBSE = tic;

startD = tic;

eqp_diff = zeros(nvbands*ncbands,1);
ev_diff = zeros(nvbands*ncbands,1);

eqp_v = eqp(nv:-1:nv-nvbands+1);
eqp_c = eqp(nvbands+1:nbands);

ev_v = ev(nv:-1:nv-nvbands+1);
ev_c = ev(nvbands+1:nbands);

for i = nv+1:nbands
  eqp_diff((i-nv-1)*nvbands+1:(i-nv)*nvbands) = eqp_c(i-nv) - eqp_v;
  ev_diff((i-nv-1)*nvbands+1:(i-nv)*nvbands) = ev_c(i-nv) - ev_v;
end

eqp_diff = eqp_diff/ry2ev;
ev_diff = ev_diff/ry2ev;

D = spdiags(eqp_diff,0,nvbands*ncbands,nvbands*ncbands);

omegalda = spdiags(ev_diff,0,nvbands*ncbands,nvbands*ncbands);

timeforD = toc(startD)

startVA = tic;

options.bse='vc';

psivcr =  prod_states_bse((psir(:,nv-nvbands+1:nv)),psir(:,nv+1:nbands),options);
psivc = F * psivcr;

Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
%Dcoul(1,1) = GWinfo.coulG0;

%GWinfo.coulG(:,4).* psivc/(8*pi);
VA = psivc' * Dcoul * psivc/(8.0*pi);
%VA_bead = psivc'*psivc
%VB = psivc' * Dcoul * conj(psivc)/(8.0*pi)

timeforVA = toc(startVA)

%chi0 = zeros(ng);
%
%startchi0 = tic;
%
%psivcr = prod_states_gw(conj(psir(:,nv-nvbands+1:nv)),psir(:,nv+1:nbands),ev,nv);
%Mg = F*psivcr;
%chi0 =  4.0 * Mg*Mg'/vol;
%
%timeforchi0 = toc(startchi0)
%
%startW = tic;
%
%I = eye(ng);
%%chi0(logical(eye(ng)))=diag(chi0).*GWinfo.coulG(:,4);
%%chi0_coul=chi0;
%%inveps = inv(I + chi0_coul)
%%chi0.*GWinfo.coulG(:,4)
%Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
%%Dcoul*chi0
%inveps = inv(I + Dcoul*chi0);
%%inveps-I
%%inveps = epsilon\I;
%%diag(inveps)
%
%W = inveps * Dcoul/8.0/pi;
%
%timeforW = toc(startW)

startWA = tic;

options.bse = 'cc';
psicc = F * prod_states_bse((psir(:,nv+1:nbands)),psir(:,nv+1:nbands),options);
options.bse = 'vv';
psivv = F * prod_states_bse((psir(:,nv-nvbands+1:nv)),psir(:,nv-nvbands+1:nv),options);
WA_tmp = psicc' * W * psivv;

WA_head_tmp = psicc(1,:)'*psivv(1,:);

WA_body = zeros(nvbands*ncbands,nvbands*ncbands);
WA_head = zeros(nvbands*ncbands,nvbands*ncbands);
WA_wing = zeros(nvbands*ncbands,nvbands*ncbands);
for i = 1:nvbands
  for j = 1:ncbands
    body_tmp = WA_tmp((j-1)*ncbands+1:j*ncbands,(i-1)*nvbands+1:i*nvbands);
    head_tmp = WA_head_tmp((j-1)*ncbands+1:j*ncbands,(i-1)*nvbands+1:i*nvbands);
    WA_body(:,i+(j-1)*nvbands) = reshape(body_tmp',nvbands*ncbands,1);
    WA_head(:,i+(j-1)*nvbands) = reshape(head_tmp',nvbands*ncbands,1);
  end
end


%WA_head
%WA_wing
%WA_body


timeforWA = toc(startWA)

%A = D+2.0*VA-WA_body

%B = zeros(nvbands*ncbands,nvbands*ncbands);

%Hbse = hbse_mat_exchange*VA+hbse_mat_body*WA_body;
%Hbse(1,2:end) = 0;
%Hbse(2:end,1) = 0;
%Hbse(1,1) = hbse_mat_head;

%omega_lda
%A = omegalda+Hbse;

%A1 = D + Hbse;

starthbse = tic;

%%calculate dipole (S0) by using momentum operator
%pol = [0,0,1];
%fac= pol*(bdot*G');
%lpol = sqrt(pol*bdot*pol');
%options.bse = 'vc';
%psivc_1 = prod_states_bse(conj(Z(:,nv-nvbands+1:nv)),Z(:,nv+1:nbands),options)/vol;
%tmp=psivc_1'.*fac;
%S = 2.0*sum(tmp,2)/lpol;
%S0 = inv(omegalda*ry2ev)*S;

%construct Hamiltonian of BSE

hbse_a =D+ hbse_mat_head*WA_head+hbse_mat_wing*WA_wing+hbse_mat_body*WA_body+hbse_mat_exchange*VA;

timeforhbse = toc(starthbse)

startdiag = tic;

[V,D,W] = eig(hbse_a);
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
Ws = W(:,ind);
%Ws'*hbse_a
%V_1 = fliplr(Vs);
%D_1 = rot90(Ds,2);
%d=real(d);
%cs = abs(Vs*conj(S0)).^2;
%hbse_a*V_1
%V(nvbands*ncbands:-1:1)
%V(end:-1:1)*conj(S0)
%V
%D

timefordiag = toc(startdiag)

%A = hbse_a;
%B = zeros(nvbands*ncbands,nvbands*ncbands);
%A1 = 0;

%startnoeh = tic;
%
%%Calculate the absorption without electron-hole interaction,
%%using eq. (27) of Rohlfing & Louie PRB 62 4927 (2000).
%emax = max(eqp_diff);
%emin = 0.0;
%iemax = ceil(emax*ry2ev+10.0*eta_abs);
%nwstep = (0:iemax/domega)/(iemax/domega);
%nwstep = repmat(nwstep,nvbands*ncbands,1);
%omega = emin + (iemax-emin)*nwstep;
%fac = omega/ry2ev-eqp_diff;
%fac_tmp = fac.^2+(eta_abs/ry2ev)^2;
%fac1 = -fac/pi./fac_tmp;
%sum1_tmp = abs(S0).^2.*fac1;
%sum1 = sum(sum1_tmp);
%fac2 = exp(-fac.^2/(2.0*(eta_abs/ry2ev)^2))/(sqrt(2.0*pi)*eta_abs/ry2ev);
%sum2_tmp = abs(S0).^2.*fac2;
%sum2 = sum(sum2_tmp);
%dos = sum(fac2);
%pref = 16*pi^2/vol;
%eps1 = 1.0 + pref*sum1;
%eps2 = pref*sum2;
%dos = dos/(ry2ev*nvbands*ncbands);
%fac = -omega/ry2ev-eqp_diff;
%fac_tmp = fac.^2+(eta_abs/ry2ev)^2;
%fac3 = -fac/pi./fac_tmp;
%sum3_tmp = abs(S0).^2.*fac3;
%sum3 = sum(sum3_tmp);
%fac4 = -exp(-fac.^2/(2.0*(eta_abs/ry2ev)^2))/(sqrt(2.0*pi)*eta_abs/ry2ev);
%sum4_tmp = abs(S0).^2.*fac4;
%sum4 = sum(sum4_tmp);
%eps3 = eps1 + pref*sum3;
%eps4 = eps2 + pref*sum4;
%abs_noeh = [omega(1,:)',eps4',eps3',dos'];
%save abs_noeh.dat -ascii abs_noeh;
%
%timefornoeh = toc(startnoeh)
%
%starteh = tic;
%
%d_tmp = real(d)*ry2ev;
%cs = abs(Vs.'*conj(S0)).^2;
%sum1 = cs'*d_tmp;
%sum1 = sum1*pref*ry2ev;
%sum1 = sum1/(0.5*pi*plasma_omega^2*ry2ev^2);
%iemax = ceil(max(d_tmp)+10.0*eta_abs);
%nwstep = (0:iemax/domega)/(iemax/domega);
%nwstep = repmat(nwstep,nvbands*ncbands,1);
%omega = emin + (iemax-emin)*nwstep;
%%d_tmp;
%fac = omega-d_tmp;
%fac_tmp = fac.^2+eta_abs^2;
%fac1 = -fac/pi./fac_tmp;
%sum1_tmp = cs.*fac1*ry2ev;
%sum1 = sum(sum1_tmp);
%fac2 = exp(-fac.^2/(2.0*(eta_abs)^2))/(sqrt(2.0*pi)*eta_abs);
%sum2_tmp = cs.*fac2*ry2ev;
%sum2 = sum(sum2_tmp);
%sum0 = sum(fac2);
%eps1 = 1.0 + pref*sum1;
%eps2 = pref*sum2;
%dos = sum0/(nvbands*ncbands);
%fac = -omega-d_tmp;
%fac_tmp = fac.^2+eta_abs^2;
%fac3 = -fac/pi./fac_tmp;
%sum3_tmp = cs.*fac3*ry2ev;
%sum3 = sum(sum3_tmp);
%fac4 = -exp(-fac.^2/(2.0*(eta_abs)^2))/(sqrt(2.0*pi)*eta_abs);
%sum4_tmp = cs.*fac4*ry2ev;
%sum4 = sum(sum4_tmp);
%eps3 = eps1+pref*sum3;
%eps4 = eps2+pref*sum4;
%abs_eh = [omega(1,:)',eps4',eps3',dos'];
%save abs_eh.dat -ascii abs_eh;
%
%timeforeh = toc(starteh)

end

timeforBSE = toc(startBSE)
