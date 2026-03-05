% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/30
%
function info = isdf_validation(what, type, dbroot, GWinfo, config)
  %
  cleanup = QPlog_push('ISDF_validation');
  % Test the density.
  
  % ---- Work space (match Fortran sizes) ----
  rho_error    = nan(1,1);
  hf_error  = nan(1,2);
  coll_error = nan(1,2);

% ------------------------------------------------------------------
% Conditional calls
  what = strtrim(what);
  if contains(what, 'rho')
    rho_error = isdf_rho_validate(type, dbroot, GWinfo);
  end
  if contains(what, 'hf')
    hf_error = isdf_hf_validate(type, dbroot, GWinfo, config);
  end
  if contains(what, 'coll')
    coll_error = isdf_coll_validate();
  end
  
% =================================================================== 
% Print test results, and warning if error is large
% ------------------------------------------------------------------
  if contains(what, 'rho')
    msg = sprintf('Density error: %.3e\n', rho_error(1));
    QPlog(msg);
    if rho_error(2) > 1e-1
      warning("Density error is large: %.3e", rho_error(1));
    end
  end

% ------------------------------------------------------------------
% Print test results, and warning if error is large
  if contains(what, 'hf')
    msg = sprintf('Hartree--Fock energy band max difference (abs, ev) = %.3e\n', hf_error(1));
    QPlog(msg, 1);
    msg = sprintf('Hartree--Fock energy band max difference (rel) = %.3e\n', hf_error(2));
    QPlog(msg, 1);
    if hf_error(1) > 0.01
      msg = sprintf(['Hartree--Fock energy band max difference %.3e~ev is too large.\n', ...
                     'Try increasing ISDF.isdf_ratio'], hf_error(1));
      warning(msg)
    end
  end

  % ---- coll ----
  if contains(what, 'coll')
    msg = sprintf('COLL difference (abs, ev) = %.3e\n', coll_error(1));
    QPlog(msg, 1);
    msg = sprintf('COLL difference (rel) = %.3e\n', coll_error(2));
    QPlog(msg, 1);
    if coll_error(2) > 0.03
      msg = sprintf(['COLL relative difference = %.3e is too large.\n', ...
                     'Try increasing ISDF.isdf_ratio'], coll_error(2));
      warning(msg);
    end
  end
end %func main
