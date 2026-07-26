function out = prod(Psi, psi, Phi, phi)
  persistent has_mex
  if isempty(has_mex)
    has_mex = ~isempty(which('isdftest.prod_mex'));
  end
  if has_mex && isa(Psi, 'double') && isa(psi, 'double') ...
      && isa(Phi, 'double') && isa(phi, 'double')
    try
      out = isdftest.prod_mex(Psi, psi, Phi, phi);
      return;
    catch ME
      if contains(ME.message, 'single only', 'IgnoreCase', true) ...
          || contains(ME.identifier, 'prod_mex:type')
        has_mex = false;
      else
        rethrow(ME);
      end
    end
  end
  out = conj(Psi * psi') .* (Phi * phi');
end
