%DEMO_ADAPTIVE_WEIGHT  Smoke demo: adaptive_weight init/get on Si relay stage.
%
% Prerequisites (from test_Si):
%   1) Run input_driver('./test') once so SAVE/config.mat and test_relay_stage.mat exist.
%   2) &ISDF isisdf = 1 in ./test
%
% Usage (from test_Si folder):
%   demo_adaptive_weight

here = fileparts(mfilename('fullpath'));
cd(here);
cd ../../;
QPstartup;
cd(here);

load('SAVE/config.mat', 'config');
if ~config.ISDF.isisdf
  error('demo_adaptive_weight:isisdf', 'Set isisdf in test input (./test) for this demo.');
end

stagePath = 'test_relay_stage.mat';
if ~exist(stagePath, 'file')
  error('demo_adaptive_weight:stage', 'Missing %s — run input_driver first.', stagePath);
end

relay.stage_from_db(stagePath);
relay.restore();

isdf.driver([], config);

L = isdf.manager('list');
id_nn = [];
for k = 1:numel(L)
  if L(k).assigned && strcmp(char(L(k).desc), 'nn')
    id_nn = L(k).id;
    break
  end
end
if isempty(id_nn)
  error('demo_adaptive_weight:id', 'No assigned ISDF slot with desc ''nn'' after isdf.driver.');
end

isdf_exclude_point(id_nn);
isdf_data = isdf.get(id_nn);
Nisdf = double(isdf_data.nisdf);
[nrange1, nrange2] = isdf_get_nrange(id_nn);
Psixga = isdf_data.coeff_seper(:, nrange1, 1, 1);
Phixga = isdf_data.coeff_seper(:, nrange2, 1, 1);

isdf_schur_update('clear');
adaptive_weight('clear');
isdf_schur_update('init', Nisdf, Psixga, Phixga);

adaptive_weight('init', id_nn);
w = adaptive_weight('get');

fprintf('demo_adaptive_weight: ISDF id = %d, desc = nn, Nisdf = %d, Nc = %d\n', ...
  id_nn, Nisdf, numel(w));
fprintf('  weight: min=%.6g  max=%.6g  sum=%.6g  mean=%.6g\n', ...
  min(w), max(w), sum(w), mean(w));

fprintf('demo_adaptive_weight: OK\n');
