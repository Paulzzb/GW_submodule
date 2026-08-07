% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/06 ZZ

function E = qp_driver(input_dir)
%QP_DRIVER  Restore SAVE services and run package QP path.
%
%   E = qp_driver(input_dir)
%
%   input_dir — SAVE directory (config.mat, data.mat, relay_stage.mat).
%   Reloads relay into memory (MATLAB-safe restart), then qp.launcher.

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

  stage_path = fullfile(input_dir, def.stage);
  if exist(stage_path, 'file') ~= 2
    error('qp_driver:stage', ...
      'Missing %s. Run input_driver first.', stage_path);
  end
  fprintf('qp_driver: loading relay stage from %s\n', stage_path);
  relay.stage_from_db(stage_path);
  relay.restore();
  if isfield(config, 'ISDF') && config.ISDF.isisdf
    fData = fullfile(input_dir, def.data);
    if exist(fData, 'file') ~= 2
      error('qp_driver:data', 'Missing %s (needed to rebuild ISDF pool).', fData);
    end
    data = load(fData, 'data').data;
    isdf.driver(data, config);
  end

  cleanup = output.push('qp_driver');
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

  output.msg('v0s', '%s', 'QP driver started (qp.launcher)');
  t0 = tic;
  E = qp.launcher(config);
  output.msg('v0s', '%s', sprintf( ...
    'quasiparticle calculation finished. total time: %.2f seconds.', toc(t0)));
end
