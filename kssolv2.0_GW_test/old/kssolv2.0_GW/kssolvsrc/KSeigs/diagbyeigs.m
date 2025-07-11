function [X,lambda] = diagbyeigs(H,neig,tol,maxitr)
%
%   [X,lambda] = diagbyeigs(H,neig,tol,maxitr) Sovle the linear 
%   eigenvalue problem HX = X*Lambda approximate by the MATLAB's eigs 
%   function
%
% Input:
%   H      the Hamiltonian object to be diagonalized
%   neig   number of occupied states (desired eigenvalues)
%   tol    convergence tolerance
%   maxitr maximum number of restarts allowed in eigs
%
n1 = H.n1;
n2 = H.n2;
n3 = H.n3;
idxnz = H.idxnz;
nnz = length(idxnz);

eigsopts.isreal = false;
eigsopts.maxit  = maxitr;
eigsopts.tol    = tol;
[V,D,flag]=eigs(@(x)eigsmult(H,x),nnz,neig,'SR',eigsopts);
d = real(diag(D));
[sd,id]=sort(d);
V = V(:,id);
if (flag~=0)
   fprintf('Convergence not reached in eigs!, pause...\');
   pause;
end
X = Wavefun(V,n1,n2,n3,idxnz);
lambda = sd;

%-----------------------------------------------------
function y = eigsmult(H,x)

n1 = H.n1;
n2 = H.n2;
n3 = H.n3;
idxnz = H.idxnz;
X = Wavefun(x,n1,n2,n3,idxnz);
Y = H*X;
y = Y.psi;
