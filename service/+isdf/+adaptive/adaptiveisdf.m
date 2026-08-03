% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/13 ZZ

function idnew = adaptiveisdf(id, cfg_isdf)
%ADAPTIVEISDF  Default adaptive path (+adaptive_double). See isdf.driver for precision routing.
  if nargin < 2
    cfg_isdf = struct();
  end
  idnew = isdf.adaptive_double.adaptiveisdf(id, cfg_isdf);
end
