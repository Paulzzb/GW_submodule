% 
% License-Identifier: BSD-3-Clause
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/30 ZZ
function validate_required_params(config)

% Validate CONTROL.groundstate_dir
if ~isfield(config, 'CONTROL') || ~isfield(config.CONTROL, 'groundstate_dir')
  msg = 'Field "CONTROL.groundstate_dir" is missing in input file.';
  output.err(msg);
end

dir_path = config.CONTROL.groundstate_dir;
if ~ischar(dir_path)
  msg = 'Field "CONTROL.groundstate_dir" must be a character vector or string.';
  output.err(msg);
end

if ~exist(dir_path, 'dir')
  msg = sprintf('Specified directory in field %s does not exist: %s', ...
        'CONTROL.groundstate_dir', dir_path);
  output.msg('v1s', '%s', msg);
  msg = sprintf('Create a directory %s in %s', ...
        'CONTROL.groundstate_dir', dir_path);
  output.msg('v1s', '%s', msg);
end

end
