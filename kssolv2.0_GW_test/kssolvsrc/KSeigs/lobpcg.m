function [X,lambda,LVEC,RVEC] = lobpcg(H, X0, prec, tol, maxit, verbose)
%
% Usage: [X,lambda,LVEC,RVEC] = lobpcg(H, X0, prec, tol, maxit, verbose);
%
% Purpose:
%    Compute the smallest eigenvalue and eigenvector of the Kohn-Sham
%    Hamiltonian using Knyazav's locally optimal preconditioned conjugate
%    gradient (LOBPCG) algorithm.
%
%    Incremental deflation against converged eigenvectors
%
% Input:
%    H     --- a Hamiltonian object
%    X0    --- the initial guess of the wave function in Fourier space
%              (Wavefun object)
%    prec  --- preconditioner
%    tol   --- tolarence on the relative residual.
%    maxit --- maximum number of iterations allowed.
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
locktol = min(1e-8,tol);
%
% get size info
%
ncol = ncols(X0); % number of wavefunctions;
if (ncol <=0) 
   fprintf('lobpcg requires at least one initial wave function!\n');
   return;
end
%
% orthonormalize the initial wave functions.
%
[X,~]=qr(X0,0); 
clear X0;
HX = H*X;
%
% some initialization merely to allocate space for W, HW, V, HV 
%
W = X;
HW = X;
V = genX0(H,ncol*3);
HV= V;
%
nconv = 0;
rcondG = 1;
%
iter = 1;
resnrm = ones(ncol,1);  % initialize residual norm
%
% --- M A I N    L O O P ---
%  
while (iter <= maxit && nconv < ncol)
  % Rayleigh quotient (approximate eigenvalue, obj func)
  S = X'*HX; S = (S+S')/2;
  [Q,D] = eig(S);
  lambda = diag(D);
  X = X*Q; HX = HX*Q;
  if ( mod(iter,10)==0 )
     HX = H*X;
  end
  V(:,1:ncol) = X;
  HV(:,1:ncol) = HX;
  R = HX - X*D;
  if (verbose == 1)
     fprintf('LOBPCG iter = %3d\n', iter);
  end 
  %
  % Check for convergence
  %
  nconv = 0;
  for j = 1:ncol 
     resnrm(j) = norm(R(:,j))/abs(lambda(j));
     %resnrm(j) = norm(R(:,j));
     if (verbose == 1)
        fprintf('eigval(%2d) = %11.3e, resnrm = %11.3e\n', j, lambda(j), resnrm(j));
     end
     if (resnrm(j) < tol) 
        nconv = nconv + 1;
     end
  end
  %
  % lock Ritz vectors that satisfy a more stringent convergence tolerance
  %
  ilock = find(abs(resnrm) <= locktol);
  iact  = find(abs(resnrm) > locktol);
  %
  nact  = length(iact);
  if (verbose == 1)
     fprintf('nconv = %d, nlock = %d\n', nconv, length(ilock));
  end;
  %
  % test for convergence
  %
  if (nconv >= ncol & iter > 1) 
     break; 
  end 
  % 
  LVEC(iter,1:ncol) = lambda.';
  RVEC(iter,1:ncol) = abs(resnrm)';
  %
  % apply the preconditioner prec
  %
  if (~isempty(prec))
     W(:,iact) = bsxfun(@times, prec, R(:,iact));
  else
     W(:,iact) = R(:,iact);
  end;
  %
  set2 = ncol+1:ncol+nact;
  W(:,iact) = W(:,iact) - X*(X'*W(:,iact));
  %[W(:,iact),~] = qr(W(:,iact),0);
  V(:,set2)  = W(:,iact);
  HW(:,iact) = H*W(:,iact); 
  HV(:,set2) = HW(:,iact);
  seta = 1:ncol+nact;
  if (iter > 1)
     set3 = ncol+nact+1:ncol+2*nact;
     V(:,set3)  = P(:,iact);
     HV(:,set3) = HP(:,iact);
     seta = 1:ncol+2*nact; 
  end;  
  %
  %  Solving Rayleigh-Ritz
  %
  T = V(:,seta)'*HV(:,seta); T = (T+T')/2;
  G = V(:,seta)'*V(:,seta); G = (G+G')/2';
  if (verbose == 1)
     fprintf('rcondG = %11.3e\n',rcond(G));
  end;
  [S,D] = eig(T,G,'chol');
  %
  U = S(:,1:ncol);
  X = V(:,seta)*U;
  HX = HV(:,seta)*U; 
  if (iter > 1)
     %
     % update the search direction
     %
     setp = ncol+1:ncol+2*nact;
     P  = V(:,setp)*U(setp,:);
     HP = HV(:,setp)*U(setp,:);
  else
     P  = W;
     HP = HW;
  end;
  %
  iter = iter + 1;
  if (verbose == 1)
     fprintf('\n'); 
  end
end
%T = X'*HX; T = (T+T')/2;
%[Q,D] = eig(T);
%X = X*Q;
%pause
