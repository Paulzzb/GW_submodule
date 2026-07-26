function init_from_config(config)
%INIT_FROM_CONFIG  Snapshot config.ISDF debug flags for the whole MATLAB session path.
%
%   Call once from service_driver (or any entry that has the authoritative GW
%   config). After this, use isdftest.debug.on(tag) and isdftest.debug.react(...) without
%   passing config into ISDF routines.
%
%   See ISDF_debug.md. Use isdftest.debug.clear with service_reset_persistent or
%   before loading a different config.

  if nargin < 1
    config = struct();
  end
  isdftest.debug.cache('init', config);
end
