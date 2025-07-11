function M = ham2mat(H)
% HAM2MAT  Convert a Hamiltonian to dense matrix.
%   M = HAM2MAT(H) converts the Hamiltonian to a dense matrix.
%
%   See also Wavefun.

n1 = H.n1;
n2 = H.n2;
n3 = H.n3;
idxnz = H.idxnz;

E = eye(numel(idxnz));
WE = Wavefun(E,n1,n2,n3,idxnz);
WM = H*WE;
TM = wavefun2tns(WM);
M = reshape(TM,n1*n2*n3,[]);
M = M(idxnz,:);

end
