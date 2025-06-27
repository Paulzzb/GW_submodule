function [X,lambda,LVEC,RVEC] = ppcg(H, X0, prec, tol, maxit, verbose)
%
% Usage: [X,lambda,LVEC,RVEC] = ppcg(H, X0, prec, tol, maxit, verbose);
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
%
% get size info
%
ncol = ncols(X0); % number of wavefunctions;
if (ncol <=0) 
   fprintf('ppcg requires at least one initial wave function!\n');
   return;
end
%
% orthonormalize the initial wave functions.
%
[X,~]=qr(X0,0);
clear X0;
HX = H*X;
R  = HX-X*(X'*HX);
%
P  = []; 
HP = []; 
%
nconv = 0;
%
maxit = 10;
iter = 1;
resnrm = ones(ncol,1);  % initialize residual norm
bsize = 2;
nb = ncol/bsize;
verbose = 1;
%
% --- M A I N    L O O P ---
%  
while (iter <= maxit && nconv < ncol)
  R = HX - X*(X'*HX);
  if (verbose == 1)
     fprintf('PPCG iter = %3d\n', iter);
  end 
  %
  % apply the preconditioner prec
  %
  W = R; 
  for j = 1:ncol
     W(:,j) = prec.*R(:,j);
  end
  %
  % orth against X, W = (I-X*X')W;
  %
  W = W - X*(X'*W);
  % 
  HW = H*W;
  %
  % construct subspaces in blocks
  %
  for jb = 1:nb
     jbeg = (jb-1)*bsize+1;
     jend = jb*bsize;
     %
     % orth sub columns of W
     %
     C = W(:,jbeg:jend)'*W(:,jbeg:jend);
     C = (C+C')/2;
     Rmat = chol(C);
     Q = W(:,jbeg:jend)/Rmat; 
     HQ = HW(:,jbeg:jend)/Rmat; 
     Y = [X(:,jbeg:jend) Q];
     HY = [HX(:,jbeg:jend) HQ];
     if (iter > 1)
        Y  = [Y P];
        HY = [HY HP];
     end;
     T = Y'*(HY); T = (T+T')/2;
     G = Y'*Y; G = (G+G')/2;
     condG = cond(G);
%    if (verbose) 
%       fprintf('iter = %d, jb = %d, cond(G) = %11.3e\n', iter, jb, condG);
%    end
     [S,D] = eig(T, G,'chol');
     [~,id]= sort(real(diag(D)));
     U = S(:,id(1:bsize));
     X(:,jbeg:jend) = Y*U;
     HX(:,jbeg:jend) = HY*U; 
     if (iter > 1)
        set2 = bsize+1:2*bsize;
        set3 = 2*bsize+1:3*bsize; 
        P(:,jbeg:jend) = W(:,jbeg:jend)*U(set2,:)  + P(:,jbeg:jend)*U(set3,:);
        HP(:,jbeg:jend) = HW(:,jbeg:jend)*U(set2,:) + HP(:,jbeg:jend)*U(set3,:);
       %
     end
  end;
  if (iter <= 1)
     %
     % initialize previous search direction
     %
     P  = W;
     HP = HW;
     for jb = 1:nb
        C = P(:,jbeg:jend)'*P(:,jbeg:jend); C = (C + C')/2;
        C = (C+C')/2;
        Rmat = chol(C);
        P(:,jbeg:jend) = P(:,jbeg:jend)/Rmat;
        HP(:,jbeg:jend) = HP(:,jbeg:jend)/Rmat;
     end;
  end;
  %
  % perform Cholesky QR on X (may choose to do it every few iterations);
  %
  C = X'*X; C = (C + C')/2;
  Rmat = chol(C);
  X = X/Rmat;
  HX = HX/Rmat;
  iter = iter + 1;
  if (verbose == 1)
     fprintf('\n'); 
  end
  %
  % check accuracy only (not to be implemented in full production)
  %
  if (0) 
     S = X'*HX; S = (S+S')/2;
     [Q,D] = eig(S);
     Z = X*Q(:,id);
     HZ = HX*Q(:,id);
     RZ = HZ-Z*D;
     for j = 1:ncol
        resnrm(j) = norm(RZ(:,j));
        fprintf('ev(%d) = %15.6e, resnrm = %11.3e\n', j, D(j,j), resnrm(j));
     end;
  end;
 % pause;
end
%
% global Rayleigh-Ritz
%
S = X'*HX; S = (S+S')/2;
[Q,D] = eig(S);
[lambda,id] = sort(real(diag(D)));
X = X*Q(:,id);
HX = HX*Q(:,id);
% check eigenvalue and residuals
R = HX-X*D(id,id);
for j = 1:ncol
   resnrm(j) = norm(R(:,j));
   fprintf('ev(%d) = %15.6e, resnrm = %11.3e\n', j, lambda(j), resnrm(j));
end;
fprintf('\n');
