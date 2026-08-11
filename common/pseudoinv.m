% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/01/27 ZZ

function Adagger = pseudoinv(A, tol)
  % Pseudoinverse of a matrix, with tolerance tol
  % Compute the SVD of A, 
  %     A = \sum \sigma_i u_i v_t'.
  % Then exclude all small singular values that sigma(i)/sigma(1) < tol,
  % Then compute 
  %     A^+ = \sum \sigma_i^{-1} v_i u_t'.
  %
  if nargin < 2
    tol = 1e-12;
  end

  [U, S, V] = svd(A);
  s = diag(S);
  cutoff = 0;
  for i = 1:size(S, 1)
    if s(i) / s(1) > tol
      cutoff = i-1;
    end
  end
  Adagger = V(:, 1:cutoff) * diag(1./s(1:cutoff)) * U(:, 1:cutoff)';
  %
end % function