function out = prod(Psi, psi, Phi, phi)
  % Validate matrix dimensions to ensure out = (Psi * psi') .* (Phi * phi') is valid
  if size(Psi, 2) ~= size(psi, 2)
      error('The number of columns in Psi (%d) must match the number of elements in psi (%d).', size(Psi,2), numel(psi));
  end
  if size(Phi, 2) ~= size(phi, 2)
      error('The number of columns in Phi (%d) must match the number of elements in phi (%d).', size(Phi,2), numel(phi));
  end
  if size(Psi,1) ~= size(Phi,1)
      error('Psi and Phi must have the same number of rows.');
  end

  % for i = 1:n1
  %   tmp1 = Psi(i, :) .* psi;
  %   tmp2 = Phi(i, :) .* phi;
  % end
  Psi = conj(Psi);
  psi = conj(psi);
  out = (Psi * psi') .* (Phi * phi');
end