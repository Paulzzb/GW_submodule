% Smoke test: adaptiveisdf(id_nn) on GaN relay stage.
%
% Prerequisites:
%   1) ./test has &ISDF isisdf = 1 (see test_GaN/test).
%   2) Run input_driver('./test') once so SAVE/config.mat and test_relay_stage.mat exist.
%
% Run from test_GaN folder:
%   run_adaptiveisdf_smoke

here = fileparts(mfilename('fullpath'));
cd(here);
cd ../../;
QPstartup;
cd(here);

load('SAVE/config.mat', 'config');
if ~config.ISDF.isisdf
  error('run_adaptiveisdf_smoke:isisdf', ...
    'Set isisdf in test input (./test) for this check.');
end

stagePath = 'test_relay_stage.mat';
if ~exist(stagePath, 'file')
  error('run_adaptiveisdf_smoke:stage', ...
    'Missing %s. Run input_driver first.', stagePath);
end

relay.stage_from_db(stagePath);
relay.restore();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isdf.driver([], config);
system_data = system.get();
wf_data = wave_functions.get();
k_data = lattice.manager('k', 'get');
nb = wf_data.nb;
nkibz = k_data.nibz;

cfg = config.ISDF;
nocc_max = 0;
nspin = system_data.nspin;
for ispin = 1:nspin
  for ikibz = 1:nkibz
    f_ib = system_data.f(:, ikibz, ispin);
    idx_last = find(f_ib(:) > 1e-5, 1, 'last');
    if ~isempty(idx_last)
      nocc_max = max(nocc_max, idx_last);
    end
  end
end

% Keep in mind that,
%     C_{n1k1n2k2, ga} = \conj{coeff_seperated(n1k1, ga)} * coeff_seperated(n2k2, ga)
% coeff_seperated: size nisdf * nb * nkibz * nspin.

isdf.free();
% id_vn = isdf.isdf_add("vn");
% vn_data = isdf.get(id_vn);
% nmu_target = cfg.isdf_ratio_type2 * sqrt(double(nocc_max) * double(nb));
% vn_data.nisdf = int32(max(1, ceil(nmu_target)));
% isdf.save2mod(vn_data, id_vn);
% %
% % Coarse-grid ISDF pipeline (ignore config.ISDF.exxmethod stubs such as kmeans).
% isdf.gen_coeff(cfg, id_vn);
% isdf.print_coarse_grid_report(id_vn);
% %
% isdf.gen_tildeVq(id_vn);
% %
% isdf.isdf_validation(id_vn);



id_nn = isdf.isdf_add("nn");
nn_data = isdf.get(id_nn);
nmu_target = cfg.isdf_ratio_type3 * double(nb);
nn_data.nisdf = int32(max(1, ceil(nmu_target)));
nn_data.assigned = true;
isdf.save2mod(nn_data, id_nn);
%
isdf.gen_coeff(cfg, id_nn);
isdf.print_coarse_grid_report(id_nn);

L = isdf.manager('list');
id_nn = [];
for k = 1:numel(L)
  if L(k).assigned && strcmp(char(L(k).desc), 'nn')
    id_nn = L(k).id;
    break;
  end
end
if isempty(id_nn)
  error('run_adaptiveisdf_smoke:id', ...
    'No assigned ISDF slot with desc ''nn'' after isdf.driver.');
end

isdf_schur_update('clear');
adaptive_weight('clear');
adaptiveisdf(id_nn);

% adaptiveisdf allocates a NEW ISDF id (idnew) for the adaptive object and runs
% isdf_validation(idnew). The HF text report is isdf_validate_HF_id<idnew>.txt
% in this folder (pwd while this script runs), not id_nn.
rep = dir(fullfile(pwd, 'isdf_validate_HF_id*.txt'));
if isempty(rep)
  warning('run_adaptiveisdf_smoke:noReport', ...
    ['No isdf_validate_HF_id*.txt under pwd=%s. ', ...
     'If adaptiveisdf finished, check earlier errors in isdf_validation.'], pwd);
else
  fprintf('run_adaptiveisdf_smoke: HF report(s) under pwd=%s :\n', pwd);
  for ri = 1:numel(rep)
    fprintf('  %s\n', fullfile(pwd, rep(ri).name));
  end
end

isdf_data = isdf.get(id_nn);
fprintf('run_adaptiveisdf_smoke: coarse slot id=%d, desc=nn, Nisdf=%d\n', ...
  id_nn, double(isdf_data.nisdf));
fprintf('run_adaptiveisdf_smoke: OK\n');
