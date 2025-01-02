function [Ex] = tddft_couplings_davidson(ksinfo, mol, nvbands, ncbands, nroots, max_dav_subspace)
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
ncv     = ncbands*nvbands;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare transition pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Normalize the wavefunction in the Fourier space
% sum_G |Psi(G)|^2 = 1
for iv = 1 : nv+ncbands
  Z(:,iv) = Z(:,iv) / (norm(Z(:,iv)));
end

% Save occupied-unoccupied pairs: psi_a*psi_i
Z_r=F'*Z(:,nv-nvbands+1:nv+ncbands);
Zcv_r = zeros(nr,ncv);
for iv = 1:nvbands
  for ic = 1:ncbands
    icv = (iv-1)*ncbands+ic;
    Zcv_r(:,icv) = conj((Z_r(:,nvbands+ic).*conj(Z_r(:,iv))));
  end
end
timeforZcv = toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build e-e interaction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coulomb kernel in reciprocal space
Dcoul_g = ksinfo.coulG;
Dcoul_mat_g = spdiags(ksinfo.coulG,0,ng,ng);
Fac_J = 2.0*vol.^3/nr;

% Exchange-correlation kernel in real space
tic
Fxc = getfxc(mol,rho);
Fxc_r = reshape(Fxc,nr,1);
Fxc_mat_r = spdiags(Fxc_r,0,nr,nr);
Fac_XC = 2.0*vol.^3/nr;

timefor2E = toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagonal part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Ediff = zeros(ncv);
for iv = 1 : nvbands
  for ic = 1 : ncbands
    icv = (iv-1)*ncbands+ic;
    Ediff(icv) = ev(nv+ic) - ev(nv-nvbands+iv);
  end
end
Ediff_mat = spdiags(Ediff,0,ncv,ncv);
[Ed,id] = sort(real(Ediff));
timefordiagonal = toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% Trial vectors
Xv = zeros(ncv,max_dav_subspace);
for iroot = 1 : nroots
  Xv(id(iroot),iroot) = 1.0;
end

% Av = H*Xv
Av = zeros(nr,max_dav_subspace);

% M = Xv'*H*Xv
M = zeros(max_dav_subspace);

max_dav_iter = 99;

ns = 0;
nl = nroots;

tol  = 1e-8;
tol2 = 1e-6;
timeforinit = toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Davidson iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter = 1 : max_dav_iter
  
  tic
  
  % Response density
  rhox  = Xv(:,ns+1:ns+nl);
  rho_r = Zcv_r*rhox;
  rho_g = F*rho_r;

  timeforDensity = toc

  tic  
  
  % \delta V_scf
  JRho_g  = (Dcoul_mat_g*rho_g);
  JRho_r  = F'*JRho_g;
  XCRho_r = Fxc_mat_r*rho_r;  
  %if(ltddft<1)
  Av(:,ns+1:ns+nl) = JRho_r*Fac_J + XCRho_r*Fac_XC;
  %elseif(ltddft<2)
  %Av(:,ns+1:ns+nl) = JRho_r*Fac_J;
  %else
  %Av(:,ns+1:ns+nl) = 0;
  %end
  
  timeforDeltaV = toc
  
  tic
  
  % Xv'*H*Xv
  Mv = rho_r'*Av(:,1:ns+nl) + rhox'*(Ediff_mat*Xv(:,1:ns+nl));  
  %Mv = rho_r'*Av(:,1:ns+nl); 
  
  timeforMv = toc
   
  for i = 1 : nl
    for j = 1 : ns+nl
      M(i+ns,j) = Mv(i,j);
      M(j,i+ns) = conj(Mv(i,j));
    end
  end  
  M2 = M(1:ns+nl,1:ns+nl);
 
  tic  
 
  % Diagonalization in subspace
  [Cx,E] = eig(M2);  
  
  timeforEig = toc

  diagE = diag(E);
  
  tic
  
  [omega,id2] = sort(real(diagE));
  
  timeforSort = toc

  Xs = Cx(:,id2);
  %omega(1:10)
  %fprintf('energy\n');
  omega_mat = spdiags(omega(1:nroots),0,nroots,nroots);
  
  tic
 
  % Build residues
  Res = Av(:,1:ns+nl)*Xs(:,1:nroots);
  Res2 = Xv(:,1:ns+nl)*Xs(:,1:nroots);
  
  timeforRes = toc
 
  tic
 
  % Project residues on CV space
  Rescv = Zcv_r'*Res + Ediff_mat*(Xv(:,1:ns+nl)*Xs(:,1:nroots)) - Res2*omega_mat;
  
  timeforRescv = toc
  
  %check = (Xv(:,1:ns+nl)*Xs(:,1:nroots))'*Rescv
  for i = 1 : nroots
    for j = 1 : ncv
      Rescv(j,i) = Rescv(j,i)/(omega(i)-Ediff(j));
    end
  end
  
  ns = ns + nl;

  for i = 1 : max_dav_subspace
    normXv(i) = Xv(:,i)'*Xv(:,i);
  end

  tic
 
  % Expand the subspace
  nl = 0;
  maxD = 0.0;
  for i = 1 : nroots
    D_Res = norm(Rescv(:,i))/ncv; 
    maxD = max(maxD,D_Res);
    if(D_Res > tol2)
      for j = 1 : ns+nl
        Rescv(:,i) = Rescv(:,i) - (Xv(:,j)'*Rescv(:,i))/normXv(j)*Xv(:,j);
      end
      normRescv  = norm(Rescv(:,i)); 
      if(normRescv > tol) 
        nl = nl + 1;
        Xv(:,ns+nl) = Rescv(:,i)/sqrt(Rescv(:,i)'*Rescv(:,i));
        normXv(ns+nl) = Xv(:,ns+nl)'*Xv(:,ns+nl);
      end
    end
  end

  timeforExpandSubspace = toc

  %fprintf('XX %3d\n',nl);
  %Xt = Xv(:,1:ns+nl);
  %XX = Xt'*Xt
  if(maxD < tol2) 
    fprintf('Converged\n');
    Ex = omega(1:nroots);
    break; 
  end
  
  if(ns+nl > max_dav_subspace)
    Xv(:,1:nroots) = Xv(:,1:ns)*Xs(:,1:nroots);
    Av(:,1:nroots) = Av(:,1:ns)*Xs(:,1:nroots);
    ns = nroots; 
  end
  
  if(nl == 0)
    fprintf('Expand subspace failed\n');
    stop;
  end

  fprintf('Iter: %3d, Nleft: %3d, Error: %10.6e\n',iter,nl,maxD);
  fprintf('time %10.6e for iter %3d\n',toc,iter)
  iter = iter + 1;  
end

if(iter > max_dav_iter)
  fprintf('maxium davidson iterations reached\n');
end

Ex = 2.0*Ex;
