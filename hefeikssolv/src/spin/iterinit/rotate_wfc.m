function [X_out, ev] = rotate_wfc(mol, H, X)
% find the best wavefunctions from the subspace which
% is spaned by the initial random wavefunctions
HX = H*X;
hc = X'*HX;
sc = X'*X;
eigsopts.isreal = false;
eigsopts.maxit  = 300;
eigsopts.tol    = 1e-16;
[V,D,flag] = eigs(hc,sc,ncols(X),'SR',eigsopts);
d = real(diag(D));
[sd,id]=sort(d);
ev = sd;
V = X.psi*V(:,id);

if (flag~=0)
   fprintf('Convergence not reached in eigs!, pause...\');
   pause;
end
X_out = Wavefun(V,mol.n1,mol.n2,mol.n3,H.idxnz,X.ispin);
end

