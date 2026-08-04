% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function E = qp_driver(input_dir)
%QP_DRIVER  Load SAVE config and run package COHSEX QP path.
%
%   E = qp_driver(input_dir)
%
%   input_dir — directory containing config.mat (usually CASE/SAVE after
%   input_driver). Delegates to qp_cohsex.

  if nargin < 1 || isempty(input_dir)
    error('qp_driver:input_dir', 'input_dir is required (e.g. ''./SAVE'').');
  end
  input_dir = char(string(input_dir));

  def = filename_map();
  fName = fullfile(input_dir, def.config);
  if exist(fName, 'file') ~= 2
    error('qp_driver:config', 'Missing %s. Run input_driver first.', fName);
  end
  TEMP = load(fName, 'config');
  config = TEMP.config;

  cleanup = output.push('qp_driver'); %#ok<NASGU>
  output.showtag(true);
  if isfield(config, 'CONTROL') && isfield(config.CONTROL, 'log_level')
    output.verbose(config.CONTROL.log_level);
  else
    output.verbose(1);
  end
  if isfield(config, 'CONTROL') && isfield(config.CONTROL, 'log_file') ...
      && ~isempty(config.CONTROL.log_file)
    output.set_logfile(config.CONTROL.log_file);
  end

  output.msg('v0s', '%s', 'QP driver started (qp_cohsex)');
  t0 = tic;
  E = qp_cohsex(config);
  output.msg('v0s', '%s', sprintf( ...
    'quasiparticle calculation finished. total time: %.2f seconds.', toc(t0)));
end
