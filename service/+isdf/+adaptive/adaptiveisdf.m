function idnew = adaptiveisdf(id, cfg_isdf)
%ADAPTIVEISDF  Default adaptive path (+adaptive_double). See isdf.driver for precision routing.
  if nargin < 2
    cfg_isdf = struct();
  end
  idnew = isdf.adaptive_double.adaptiveisdf(id, cfg_isdf);
end
