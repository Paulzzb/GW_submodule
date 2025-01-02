function [X,lambda] = davidson2_gpu(H, X0, prec, tol, maxit, verbose)
%
% Usage: [X,lambda] = davidson(H, X0, prec, tol, maxit, verbose);
%
% Purpose:
%    Compute the smallest eigenvalue and eigenvector of the Kohn-Sham
%    Hamiltonian using Davidson iterative algorithm 
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
%

X      = [];
lambda = [];

if ~isa(H, 'function_handle'), H = @(x)H*x; end

%
% get size info
%
if isa(X0,'Wavefun')
   ncol = ncols(X0); % number of wavefunctions;
else
   ncol = size(X0,2);
end
if (ncol <=0) 
   fprintf('davidson requires at least one initial wave function!\n');
   return;
end
%
% orthonormalize the initial wave functions.
%
X0.psi = gpuArray(X0.psi);
[X,~]=qr(X0,0);
%clear X0;
HX = H(X);
%
nconv = 0;
%
iter = 1;
resnrm = gpuArray.ones(ncol,1);  % initialize residual norm
%
% --- M A I N    L O O P ---
%  
etol = gpuArray.zeros(ncol,1);
etol(1:ncol) = max((tol*5.0),1e-5);
mtype = find(X.occ==1);
etol(mtype) = tol;
%iact = 1:ncol;
%nact = ncol;
eold = gpuArray.zeros(ncol,1);
while (iter <= maxit && nconv < ncol)
  % Rayleigh quotient (approximate eigenvalue, obj func)
  if (iter==1)
    S = X'*HX;
  end
  lambda = eig(S);
  lambda = sort(real(lambda)); 

  R = HX - X*S;
  if (verbose == 1)
     fprintf('DAVIDSON iter = %3d\n', iter);
  end 
  %
  % Check for convergence
  %
  for j = 1:ncol 
     resnrm(j) = norm(R(:,j));
     %if (verbose == 1)
     %   fprintf('eigval(%2d) = %11.3e, resnrm = %11.3e\n', j, lambda(j), resnrm(j));
     %end
  end
  iconv = find(abs(lambda(1:ncol)-eold(1:ncol))<=etol);
  iact  = find(abs(lambda(1:ncol)-eold(1:ncol))> etol);
  %iact  = find(abs(resnrm) > etol);
  % 
  nconv = length(iconv);
  nact  = length(iact);
  fprintf('nact = %d, nconv = %d\n', nact,nconv);
  %
  % test for convergence
  %
  if (nconv >= ncol) 
     break; 
  end 
  eold  = lambda;
  %
  % apply the preconditioner prec
  %
  W = bsxfun(@times,prec,R(:,iact));
  W = W - X*(X'*W);
  for j = 1:nact
     W(:,j) = W(:,j)/norm(W(:,j));
  end
  
  HW = H(W);
  Q  = [X W];
  HQ = [HX HW];

  %T = Q'*(HQ)
  %G = Q'*Q
  T1 = diag(real(lambda(1:ncol)));
  T2 = W'*(HX);
  T3 = W'*(HW);
  T  = [T1 T2'; T2 T3];
  T  = (T+T')/2;
  G1 = gpuArray.eye(ncol);
  G2 = W'*X;
  G3 = W'*W;
  G  = [G1 G2'; G2 G3];
  G = (G+G')/2;
  condG = cond(G);
  if (condG > 1.0e12) 
     fprintf('cond(G) = %11.3e\n', condG);
     break;
  end

 % tic
 % T = gather(T);
 % G = gather(G);
 % [S,D] = eig(T, G,'chol');
 % T = gpuArray(T);
 % G = gpuArray(G);
 % S = gpuArray(S);
 % D = gpuArray(D);
  invG = inv(G);
  Eig = invG*T;
  [S,D] = eig(Eig);

 % t_eig = toc

  [~,id]= sort(real(diag(D)));
  U = S(:,id(1:ncol));
  X = Q*U;
  HX = HQ*U; 
  S = U'*(T*U);
  %
  iter = iter + 1;
  if (verbose == 1)
     fprintf('\n'); 
  end
  %pause
end
S = X'*HX; S = (S+S')/2;
[Q,D] = eig(S);
[lambda,id] = sort(real(diag(D)));
X = X*Q(:,id);

X.psi  = gather(X.psi);
lambda = gather(lambda);
end
