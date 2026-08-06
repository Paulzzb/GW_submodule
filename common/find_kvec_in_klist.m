% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/19 ZZ

function [ind_k] = find_kvec_in_klist(kpt, kpt_set, TOL)
  % Return the index of kpt in kpt_set
  %       -1 if not found
  % Default TOL = 1e-9
  if nargin < 3
    TOL = 1e-5;
  end
  
  ind_k = -1;
  nbz = length(kpt_set(:)) / 3;

  
  for ii=1:nbz
    kpt_diff = kpt - kpt_set(ii,:);
    v = kpt_diff - round(kpt_diff);
    if norm(v) <= TOL 
      ind_k = ii;
      return;
    end
  end

  warning("find_k_indx: kpt not found in kpt_set");
end % EOF 




