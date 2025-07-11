function [X,lambda,LVEC,RVEC] = rmmdiis(H, X0, prec, tol, maxdiis, verbose)
%
% Usage: [X,lambda,LVEC,RVEC] = rmmdiis(H, X0, prec, tol, maxdiis, verbose);
%
% Purpose:
%    Compute the smallest eigenvalue and eigenvector of the Kohn-Sham
%    Hamiltonian using the residual minimization method (RMM) combined
%    with direct inversion of iterative subspace (DIIS)
%
%    Incremental deflation against converged eigenvectors
%
% Input:
%    H       --- a Hamiltonian object
%    X0      --- the initial guess of the wave function in Fourier space
%               (Wavefun object)
%    prec    --- preconditioner
%    tol     --- tolarence on the relative residual.
%    maxdiis --- maximum dimension of the DIIS subspace
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

X      = X0;
%
% get size info
%
ncol = ncols(X0); % number of wavefunctions;
if (ncol <=0) 
   fprintf('Davidson requires at least one initial wave function!\n');
   return;
end
lambda = zeros(ncol,1) ;
LVEC   = zeros(maxdiis+1,ncol);
RVEC   = zeros(maxdiis+1,ncol);
%
V  = genX0(H,maxdiis);  % temporary storage for DIIS vectors
W  = V;                 % needed to hold DIIS difference vectors
HV = V;                 % needed to update Hvb without multiplying H explicitly
for j = 1:ncol
   v = X0(:,j)/norm(X0(:,j));
   Hv = H*v;
   V(:,1) = v;
   HV(:,1) = Hv;
   theta = v'*Hv; 
   r = Hv - theta*v; 
   Hvb = Hv;
   vb = v;
   w1 = bsxfun(@times, prec, r);
   resnrm = norm(r)/abs(theta); 
   RVEC(1,j) = resnrm;
   iter = 1;
   fprintf('iter = %d, theta = %15.6e, resnrm = %11.3e\n', iter, theta, resnrm);
   while (iter <= maxdiis && resnrm > tol)
      if (iter > 1)
         w = bsxfun(@times, prec, r);
         W(:,iter-1) = w1-w;
         a = W(:,1:iter-1)\w1;
         a1 = 1-sum(a);
         vb = V(:,1)*a1 + V(:,2:iter)*a;
         nrmvb = norm(vb);
         vb = vb/nrmvb;
         Hvb = HV(:,1)*a1/nrmvb+HV(:,2:iter)*a/nrmvb;
         theta = vb'*Hvb;
      end

      rb = Hvb - theta*vb;
      wb = bsxfun(@times, prec, rb);
      Q = [vb wb];
      HQ = [Hvb H*wb];
      T = Q'*HQ; G = Q'*Q;
      T = (T+T')/2; G = (G+G')/2;
      [S,D]=eig(T,G,'chol');
      v = Q*S(:,1);
      theta = D(1,1);
      Hv = HQ*S(:,1);
      V(:,iter+1) = v;
      HV(:,iter+1) = Hv;
      %
      r = Hv-v*theta;
      resnrm = norm(r)/abs(theta);
      iter = iter + 1;
      RVEC(iter,j) = resnrm;
      fprintf('iter = %2d, theta = %15.6e, resnrm = %11.3e\n', iter, theta, resnrm);
      if (resnrm < tol) 
         break;
      end
   end;
   X(:,j) = v;
   lambda(j) = theta;
   fprintf('\n');
   disp('pause');
   pause
end;
%
