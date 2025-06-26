function validate_required_params(config)

% Validate CONTROL.groundstate_dir
if ~isfield(config, 'CONTROL') || ~isfield(config.CONTROL, 'groundstate_dir')
  msg = 'Field "CONTROL.groundstate_dir" is missing in input file.';
  QPerror(msg);
end

dir_path = config.CONTROL.groundstate_dir;
if ~ischar(dir_path)
  msg = 'Field "CONTROL.groundstate_dir" must be a character vector or string.';
  QPerror(msg);
end

if ~exist(dir_path, 'dir')
  msg = sprintf('Specified directory in field %s does not exist: %s', ...
        'CONTROL.groundstate_dir', dir_path);
  QPerror(msg);
end

end
