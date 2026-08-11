% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06 ZZ

function data = load_kssolv_groundstate(dirin)
% load_kssolv_groundstate - Load groundstate data from a KSSOLV-style .mat file
% User must have previously run a KSSOLV calculation and saved the results
% through save_groundstate_to_GWformat.m

filePath = fullfile(dirin, 'groundstate.mat');

if ~exist(filePath, 'file')
  msg = sprintf('groundstate.mat not found in directory: %s', dirin);
  output.err(msg);
end

tmp = load(filePath);

if ~isfield(tmp, 'groundstate')
  msg = 'groundstate.mat does not contain a variable named "groundstate".';
  output.err(msg);
end

required_fields = {'rhor', 'Vxc', 'ev', 'psig', 'sys', 'occupation', 'reciprocal_grid_info', 'nkibz', 'nspin', 'kibz', 'nspinor'};
for k = 1:length(required_fields)
  if ~isfield(tmp.groundstate, required_fields{k})
    msg = sprintf('Missing field "%s" in groundstate structure.', required_fields{k});
    output.err(msg);
  end
end

data = tmp.groundstate;

end
