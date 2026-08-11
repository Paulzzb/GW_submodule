% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06 ZZ

function init_from_config(config)
%INIT_FROM_CONFIG  Snapshot config.ISDF debug flags for the whole MATLAB session path.
%
%   Call once from service_driver (or any entry that has the authoritative GW
%   config). After this, use isdf.debug.on(tag) and isdf.debug.react(...) without
%   passing config into ISDF routines.
%
%   See ISDF_debug.md. Use isdf.debug.clear with service_reset_persistent or
%   before loading a different config.

  if nargin < 1
    config = struct();
  end
  isdf.debug.cache('init', config);
end
