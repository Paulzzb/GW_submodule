% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/23

function report = restore(stage)
% Restore staged relay data back into service managers.
% Based on relay configuration: use relay.get_relay_config to define what to restore.
% If stage is omitted, use the staged cache from appdata.

  if nargin < 1 || isempty(stage)
    if ~isappdata(0, 'GW_RELAY_STAGE')
      error('relay::restore missing staged data. Run relay.collect or relay.stage_from_db first.');
    end
    stage = getappdata(0, 'GW_RELAY_STAGE');
  end

  cfg = relay.get_relay_config();

  report = struct();
  report.ok = true;
  report.saved = {};
  report.skipped = {};
  report.errors = {};

  % Dynamically iterate over configured modules and variables
  modules = fieldnames(cfg.modules);
  for m_idx = 1:numel(modules)
    mod_name = modules{m_idx};
    mod_cfg = cfg.modules.(mod_name);
    
    var_names = fieldnames(mod_cfg);
    for v_idx = 1:numel(var_names)
      var_name = var_names{v_idx};
      var_cfg = mod_cfg.(var_name);
      
      % Only process if this has call_set (is a trackable variable)
      if isfield(var_cfg, 'call_set')
        full_name = [mod_name '.' var_name];
        if isfield(stage, mod_name) && isfield(stage.(mod_name), var_name)
          try
            var_cfg.call_set(stage.(mod_name).(var_name));
            report.saved{end + 1} = full_name;
          catch ME
            report.ok = false;
            report.errors{end + 1} = [full_name ': ' ME.message];
          end
        else
          report.skipped{end + 1} = full_name;
        end
      end
    end
  end

  setappdata(0, 'GW_RELAY_STAGE', stage);
end
