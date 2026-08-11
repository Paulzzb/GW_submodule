% Smoke test: isdf.driver([], config) after relay restore.
% Run from test_Si folder: run_isdf_driver_smoke

here = fileparts(mfilename('fullpath'));
cd(here);
cd ../../;
QPstartup;
cd(here);

load('SAVE/config.mat', 'config');
if ~config.ISDF.isisdf
  error('run_isdf_driver_smoke: set isisdf in test input for this check.');
end

stagePath = 'test_relay_stage.mat';
if ~exist(stagePath, 'file')
  error('run_isdf_driver_smoke: missing %s (run input_driver first).', stagePath);
end
relay.stage_from_db(stagePath);
relay.restore();

isdf.driver([], config);
fprintf('run_isdf_driver_smoke: OK\n');
