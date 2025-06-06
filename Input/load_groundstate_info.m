function data = load_groundstate_info(dirin, typein)
% Load groundstate information and convert to @GWinfor format

% Step 1: Read raw data
switch lower(typein)
  case 'kssolv'
    data = load_kssolv_groundstate(dirin);
  case 'qe'
    data = load_qe_groundstate(dirin);
  otherwise
    error('Unsupported groundstate type: %s', typein);
end

end


function data = load_kssolv_groundstate(dirin)
% load_kssolv_groundstate - Load groundstate data from a KSSOLV-style .mat file
% User must have previously run a KSSOLV calculation and saved the results
% through save_groundstate_to_GWformat.m

filePath = fullfile(dirin, 'groundstate.mat');

if ~exist(filePath, 'file')
    error('groundstate.mat not found in directory: %s', dirin);
end

tmp = load(filePath);

if ~isfield(tmp, 'groundstate')
    error('groundstate.mat does not contain a variable named "groundstate".');
end

required_fields = {'rhor', 'Vxc', 'ev', 'psi', 'sys', 'occupation', 'reciprocal_grid_info'};
for k = 1:length(required_fields)
    if ~isfield(tmp.groundstate, required_fields{k})
        error('Missing field "%s" in groundstate structure.', required_fields{k});
    end
end

data = tmp.groundstate;

end

function data = load_qe_groundstate(dirin)
% load_qe_groundstate - Load groundstate data from a QE pw.x output.
% User must have previously run a QE pw.x calculation.
error('load_qe_groundstate not implemented yet.');
end

