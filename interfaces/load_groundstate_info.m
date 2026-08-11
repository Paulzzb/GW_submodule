% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06 ZZ

function data = load_groundstate_info(dirin, typein, config)
% Load groundstate information and convert to @GWinfor format
%
%   data = load_groundstate_info(dirin, typein)
%   data = load_groundstate_info(dirin, typein, config)  % required for 'formal'

if nargin < 3
  config = [];
end

% Step 1: Read raw data
switch lower(typein)
  case 'kssolv'
    data = load_kssolv_groundstate(dirin);
  case 'qe'
    data = load_qe_groundstate(dirin);
  case 'formal'
    if isempty(config)
      error('load_groundstate_info:formal', ...
        'groundstate_type=''formal'' requires config (call from input_driver).');
    end
    error('load_groundstate_info:formal', ...
      'groundstate_type=''formal'' is not supported in this version.');
    data = load_formal_groundstate(dirin, config);
  otherwise
    msg = sprintf('Unsupported groundstate type: %s', typein);
    output.err(msg);
end

end

