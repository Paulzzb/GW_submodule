function [ind_k] = find_gvec_in_glist(kpt, kpt_set, fftgrid, TOL)
  % Return the index of kpt in kpt_set
  %       -1 if not found
  % Default TOL = 1e-9
  if nargin < 4
    TOL = 1e-7;
  end
  
  nbz = length(kpt_set(:)) / 3;

  flaghalf = ( abs(round(fftgrid/2) - fftgrid/2) < TOL);
  halffft = round(fftgrid / 2);
  mask = flaghalf & (kpt+TOL > halffft);
  kpt_new = kpt;
  kpt_new(mask) = kpt(mask) - fftgrid(mask);

  diff = abs(kpt_set - kpt_new);
  dist = sum(diff,2);
  ind_k = find(dist <= TOL, 1);

  if isempty(ind_k)
    ind_k = -1;
    error();
  end
  % warning("find_k_indx: kpt not found in kpt_set");
end % EOF 




