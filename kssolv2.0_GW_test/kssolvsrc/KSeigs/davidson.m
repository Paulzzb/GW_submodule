function [X,lambda,LVEC,RVEC] = davidson(H, X0, prec, tol, nsteps, maxcyc, verbose)
%
% Usage: [X,lambda,LVEC,RVEC] = davidson(H, X0, prec, tol, nsteps, maxcyc, verbose);
%
% Purpose:
%    Compute the smallest eigenvalue and eigenvector of the Kohn-Sham
%    Hamiltonian using the Davidson algorithm
%
%    Incremental deflation against converged eigenvectors
%
% Input:
%    H       --- a Hamiltonian object
%    X0      --- the initial guess of the wave function in Fourier space
%               (Wavefun object)
%    prec    --- preconditioner
%    tol     --- tolarence on the relative residual.
%    nsteps  --- number of Davidson steps before restart
%    maxcyc  --- maximum number of restart cycles
%
% Output:
%    X      --- eigenvectors corresponding to the smallest eigenvalues of H
%               (Wavefun object)
%    lambda --- approximations to the smallest eigenvalues of H
%    LVEC   --- eigenvalue history (matrix of size m by ncols, where m is the 
%               total number of iterations and ncols is the number of eigenvalues
%               computed.)
%    RVEC   --- residual norm history (matrix of size m by ncols)
%

X      = [];
lambda = [];
LVEC   = [];
RVEC   = [];
locktol = min(1e-5,tol);
%
% get size info
%
ncol = ncols(X0); % number of wavefunctions;
if (ncol <=0) 
   fprintf('Davidson requires at least one initial wave function!\n');
   return;
end
%
% orthonormalize the initial wave functions.
%
HX = H*X0;
T = X0'*HX; T = (T+T')/2;
G = X0'*X0; G = (G+G')/2;
if ( norm(G-eye(size(G)),'fro')/size(G,1) < 1e-15 )
   [S,D]=eig(T);
else
   [S,D]=eig(T,G,'chol');
end;
X = X0*S;
HX = HX*S;
clear X0;
%
nconv = 0;
%
% create some objects to hold Davidson subspace  
%
V  = genX0(H,ncol*nsteps);
HV = V;
V(:,1:ncol) = X;
HV(:,1:ncol) = HX;
HV1 = H*V(:,1:ncol);
nv = ncol;
W  = X; 

iter = 1;
resnrm = ones(ncol,1);  % initialize residual norm
%
% --- M A I N    L O O P ---
%  
while (iter <= maxcyc && nconv < ncol)
  if (verbose == 1)
     fprintf('Davidson restart cycle: %3d\n', iter);
  end 
  for itd = 1:nsteps
     if (verbose == 1)
        fprintf('   Davidson step: %3d\n', itd);
     end 
     % Rayleigh quotient (approximate eigenvalue, obj func)
     lambda = diag(D(1:ncol,1:ncol));
     R = HX - X*D(1:ncol,1:ncol);
     %
     % Check for convergence
     %
     for j = 1:ncol 
        resnrm(j) = norm(R(:,j))/abs(lambda(j));
        if (verbose == 1)
           fprintf('eigval(%2d) = %11.3e, resnrm = %11.3e\n', j, lambda(j), resnrm(j));
        end
     end
     iconv = find(abs(resnrm)<=tol);
     %
     % lock Ritz vectors that satisfy a more stringent convergence tolerance
     %
     ilock = find(abs(resnrm) <= locktol);
     iact  = find(abs(resnrm) > locktol);
     %
     % nconv is generally larger than nlock
     % 
     nconv = length(iconv);
     nlock = length(ilock);
     nact  = length(iact);
     if (verbose == 1)
        fprintf('nlock = %d, nconv = %d\n', nlock,nconv);
     end
     %
     % test for convergence
     %
     if (nconv >= ncol & itd > 1) 
        break; 
     end 
     % 
     LVEC((iter-1)+itd,1:ncol) = lambda';
     RVEC((iter-1)+itd,1:ncol) = abs(resnrm)';
     %
     % apply the preconditioner prec
     %
     if (~isempty(prec))
        W(:,iact) = bsxfun(@times, prec, R(:,iact));
     else
        W(:,iact) = R(:,iact);
     end;
     %
     b1 = 1:nv;
     b2 = nv+1:nv+nact;
     %
     W(:,iact) = W(:,iact) - V(:,b1)*(V(:,b1)'*W(:,iact));
     [V(:,b2),R] = qr(W(:,iact),0); % can be replaced by Cholesky QR
     HV(:,b2) = H*V(:,b2);
     T(b2,b2) = V(:,b2)'*HV(:,b2);
     T(b1,b2) = V(:,b1)'*HV(:,b2);
     T(b2,b1) = T(b1,b2)';
     nv = nv + nact;
     T(1:nv,1:nv) = (T(1:nv,1:nv)+T(1:nv,1:nv)')/2;
     %
     [S,D] = eig(T(1:nv,1:nv));
     X  = V(:,1:nv)*S(:,1:ncol);
     HX = HV(:,1:nv)*S(:,1:ncol);
     %
  end % for itd
  if (verbose == 1)
     fprintf('\n'); 
  end
  V(:,1:ncol) = X;
  HV(:,1:ncol) = HX;
  nv = ncol;
  T(1:nv,1:nv) = X'*HX;
  T(1:nv,1:nv) = (T(1:nv,1:nv)+T(1:nv,1:nv)')/2;
  iter = iter + 1;
  %pause
end % while
T = X'*HX; T = (T+T')/2;
[S,D] = eig(T);
[lambda,id] = sort(real(diag(D)));
X = X*S(:,id);
