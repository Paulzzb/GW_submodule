function validate_required_params(config)

% Validate CONTROL.groundstate_dir
if ~isfield(config, 'CONTROL') || ~isfield(config.CONTROL, 'groundstate_dir')
    error('[GW] Required field "CONTROL.groundstate_dir" in input file is missing.');
end

dir_path = config.CONTROL.groundstate_dir;
if ~ischar(dir_path)
    error('[GW] CONTROL.groundstate_dir must be a character vector or string.');
end

if ~exist(dir_path, 'dir')
    error('[GW] Specified directory in field %s does not exist: %s', ...
          'CONTROL.groundstate_dir', dir_path);
end

end
