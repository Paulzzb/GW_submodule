function w = ncols(BH)
% BLOCHHAM/NCOLS Numbers of columns in the Bloch Hamiltonian class
%    w = NCOLS(BH) returns the numbers of columns in the non-local pseudo
%    potential of the Bloch Hamiltonian.
%
%    See also BlochHam.

w = cellfun(@(X)size(X,2),BH.vnlmatcell);

end
