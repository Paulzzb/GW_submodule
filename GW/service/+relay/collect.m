% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/23

function stage = collect()
% Collect currently initialized persistent states from service managers.
% Based on relay configuration: use relay.get_relay_config to define what to collect.
% The collected snapshot is also cached in appdata for relay.restore.

  cfg = relay.get_relay_config();
  
  stage = struct();
  stage.meta = struct();
  stage.meta.version = cfg.version;
  stage.meta.source = 'relay.collect';
  stage.meta.collected_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
  stage.status = struct();

  % Dynamically iterate over configured modules and variables
  modules = fieldnames(cfg.modules);
  for m_idx = 1:numel(modules)
    mod_name = modules{m_idx};
    mod_cfg = cfg.modules.(mod_name);
    
    % Initialize module in stage if not present
    if ~isfield(stage, mod_name)
      stage.(mod_name) = struct();
    end
    
    var_names = fieldnames(mod_cfg);
    for v_idx = 1:numel(var_names)
      var_name = var_names{v_idx};
      var_cfg = mod_cfg.(var_name);
      
      % Only process if this has call_get (is a trackable variable)
      if isfield(var_cfg, 'call_get')
        status_key = [mod_name '_' var_name];
        try
          stage.(mod_name).(var_name) = var_cfg.call_get();
          stage.status.(status_key) = 'ok';
        catch ME
          stage.status.(status_key) = ['skip: ' ME.message];
        end
      end
    end
  end

  setappdata(0, 'GW_RELAY_STAGE', stage);
end
