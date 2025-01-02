function [X1,lambda] = davidson_cpu(H, X0, prec, tol, maxit, verbose)
%
% Usage: [X1,lambda,LVEC,RVEC] = davidson(H, X0, prec, tol, maxit, verbose);
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
%    X1     --- eigenvectors corresponding to the smallest eigenvalues of H
%               (Wavefun object)
%    lambda --- approximations to the smallest eigenvalues of H

eold   = H.eband;

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
%X0.psi = gpuArray(X0.psi);
[X,~]=qr(X0,0);
%clear X0;
HX = H(X);
%
nconv = 0;
%
iter = 1;
resnrm = ones(ncol,1);  % initialize residual norm
%
% --- M A I N    L O O P ---
%  
% initialize work matrices
Hs = X'*HX;
S  = X'*X;
Q  = X;
HQ = HX;

% Rayleigh quotient (approximate eigenvalue, obj func)
[V,D] = eig(Hs);
[lambda,id]= sort(real(diag(D)));
U = V(:,id);
R = HX - X*Hs;

nsub = ncol;
maxsub = max(4*ncol,300);

while (iter <= maxit && nconv < ncol)
  %if (verbose == 1)
  %   fprintf('DAVIDSON iter = %3d\n', iter);
  %end 
  %
  % Check for convergence
  %
  for j = 1:ncol 
     resnrm(j) = norm(R(:,j));
     %if (verbose == 1)
     %   fprintf('eigval(%2d) = %11.3e, resnrm = %11.3e\n', j, lambda(j), resnrm(j));
     %end
  end
  iconv = find(abs(lambda(1:ncol)-eold(1:ncol))<=tol);
  %iact  = find(abs(lambda(1:ncol)-eold(1:ncol))> tol);
  iact  = find(abs(resnrm) > tol);
  nconv = length(iconv);
  nact  = length(iact);
  nsub  = nsub + nact;
  %if (verbose == 1)
  %  fprintf('nact = %d, nconv = %d\n', nact,nconv);
  %end
  %
  % update eigenvalues
  %
  eold  = lambda;
  %
  % test for convergence
  %
  if (nconv >= ncol) 
     break; 
  elseif (nsub > maxsub)
     fprintf('maximum subspace!!!restart davidson\n')
     Q  = Q*U;
     HQ = HQ*U;
     Hs = U'*Hs*U;
     S  = eye(ncol);
     nsub = ncol; 
  end 
  % 
  %LVEC(iter,1:ncol) = lambda.';
  %RVEC(iter,1:ncol) = abs(resnrm)';
  %
  % apply the preconditioner prec
  %
  X = bsxfun(@times,prec,R(:,iact));
  for j = 1:nact
     X(:,j) = X(:,j)/norm(X(:,j));
  end

  G2 = X'*Q;
  G3 = X'*X;
  S  = [S G2'; G2 G3];
  condS = cond(S);
  if (condS > 1.0e12) 
     fprintf('cond(S) = %11.3e\n', condS);
     break;
  end

  HX = H(X);
  T2 = X'*(HQ);
  T3 = X'*(HX);
  Hs  = [Hs T2'; T2 T3];

  %invS = inv(S);
  %Eig = invS*Hs;
  %[V,D] = eig(Eig);
  [V,D] = eig(Hs,S,'chol');
  [lambda,id]= sort(real(diag(D)));
  U = V(:,id(1:ncol));

  Q  = [Q X];
  HQ = [HQ HX];
  R = (HQ*U) - (Q*U)*(U'*Hs*U);
  %end 
  %
  iter = iter + 1;
  %if (verbose == 1)
  %   fprintf('\n'); 
  %end
  %pause
end % while

X = Q*U;
S = U'*Hs*U; 
[V,D] = eig(S);
[lambda,id] = sort(real(diag(D)));
X1 = X*V(:,id);

%X1.psi = gather(X1.psi);
%lambda = gather(lambda);

end % function
