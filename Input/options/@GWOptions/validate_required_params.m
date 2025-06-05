function validate_required_params(config)
% Ensure all required parameters are present and valid
% See GW_input_description for a list of required parameters

% CONTROL.groundstate_dir
if ~isfield(config, 'CONTROL') || ~isfield(config.CONTROL, 'groundstate_dir')
    error('Missing required parameter: CONTROL.groundstate_dir');
end

groundstate_dir = config.CONTROL.groundstate_dir;

if isempty(groundstate_dir) || ~ischar(groundstate_dir)
    error('CONTROL.groundstate_dir must be a non-empty string.');
end

if ~exist(groundstate_dir, 'dir')
    error('The specified CONTROL.groundstate_dir does not exist: %s', groundstate_dir);
end

% CONTROL.groundstate_type must be 'kssolv' or 'qe'
if ~isfield(config.CONTROL, 'groundstate_type')
    error('Missing required parameter: CONTROL.groundstate_type');
end

valid_types = {'kssolv', 'qe'};
if ~any(strcmpi(config.CONTROL.groundstate_type, valid_types))
    error('CONTROL.groundstate_type is not supported. Currently only ''kssolv'' and ''qe'' are allowed.');
end

% Todo

end
