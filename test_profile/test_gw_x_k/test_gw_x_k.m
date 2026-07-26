cd ../../
QPstartup
cd test_profile/test_bz_samp

load ../TMP_FILES/GWinput.mat
load ../TMP_FILES/config.mat

% Some preparations for testing
filePath = 'test_relay_stage.mat';

if ~exist(filePath, 'file')
  error('test_stage: File not found: %s', filePath);
end

fprintf('\n=== Loading Relay Stage ===\n');
fprintf('Loading from: %s\n', filePath);

% Load stage from database
relay.stage_from_db(filePath);
fprintf('Stage loaded into relay cache.\n');

% Restore to all modules
fprintf('\n=== Restoring to Modules ===\n');
report = relay.restore();

% Display restore report
if report.ok
  status_str = 'OK';
else
  status_str = 'FAILED';
end

% GWgroundstate.psir = get_wavefunc_real(GWgroundstate.psig, GWgroundstate.Ggrid4psig);


Ex = gw_x_k(GWgroundstate, config);


