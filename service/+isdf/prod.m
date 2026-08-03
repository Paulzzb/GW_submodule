% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/13 ZZ

function out = prod(Psi, psi, Phi, phi)
  persistent has_mex
  if isempty(has_mex)
    has_mex = ~isempty(which('isdf.prod_mex'));
  end
  if has_mex && isa(Psi, 'double') && isa(psi, 'double') ...
      && isa(Phi, 'double') && isa(phi, 'double')
    try
      out = isdf.prod_mex(Psi, psi, Phi, phi);
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
