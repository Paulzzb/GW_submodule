function [ind_k] = find_vec_in_list(kpt, kpt_set, TOL)
  % Return the index of kpt in kpt_set
  %       -1 if not found
  % Default TOL = 1e-9
  if nargin < 3
    TOL = 1e-9;
  end

  ind_k = -1;
  nbz = length(kpt_set(:)) / 3;

  for ii=1:nbz
    tmpf=abs(kpt-kpt_set(ii,:));
    if sum(abs(tmpf)) <= TOL 
      ind_k = ii;
      return;
    end
  end

  warning("find_k_indx: kpt not found in kpt_set");
end % EOF 




