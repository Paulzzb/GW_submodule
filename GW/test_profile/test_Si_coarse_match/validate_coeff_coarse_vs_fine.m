% Validate coarse-grid coeff_seper against fine-grid wf_data.c on mapped R points.
%
% Run from this folder:
%   validate_coeff_coarse_vs_fine

here = fileparts(mfilename('fullpath'));
cd(here);
cd ../../;
QPstartup;
cd(here);

stagePath = 'test_relay_stage.mat';
if ~exist(stagePath, 'file')
  error('validate_coeff_coarse_vs_fine:stage', ...
    'Missing %s. Please ensure stage file is copied.', stagePath);
end

cfgPath = fullfile('SAVE', 'config.mat');
if ~exist(cfgPath, 'file')
  error('validate_coeff_coarse_vs_fine:config', ...
    'Missing %s. Please ensure SAVE/config.mat is copied.', cfgPath);
end

relay.stage_from_db(stagePath);
relay.restore();

load(cfgPath, 'config');
isdf.driver([], config);

% Locate current "nn" ISDF slot.
L = isdf.manager('list');
id_nn = [];
for k = 1:numel(L)
  if L(k).assigned && strcmp(char(L(k).desc), 'nn')
    id_nn = L(k).id;
    break;
  end
end
if isempty(id_nn)
  error('validate_coeff_coarse_vs_fine:id', ...
    'No assigned ISDF slot with desc ''nn'' after isdf.driver.');
end

isdf_data = isdf.get(id_nn);
wf_data = wave_functions.get();
fft_data = FFT.get();

fftgrid = int32(fft_data.fftgrid(:).');
fftgrid_c = int32(isdf_data.fftgrid_c(:).');
if any(fftgrid ~= int32([18, 18, 18]))
  error('validate_coeff_coarse_vs_fine:fineGrid', ...
    'Expected fine fftgrid=[18,18,18], got [%d,%d,%d].', ...
    fftgrid(1), fftgrid(2), fftgrid(3));
end
if any(fftgrid_c ~= int32([9, 9, 9]))
  error('validate_coeff_coarse_vs_fine:coarseGrid', ...
    'Expected coarse fftgrid_c=[9,9,9], got [%d,%d,%d].', ...
    fftgrid_c(1), fftgrid_c(2), fftgrid_c(3));
end

R_sampling_RLU = double(isdf_data.R_sampling_RLU);
Nmu = size(R_sampling_RLU, 1);
nb = size(isdf_data.coeff_seper, 2);

if size(isdf_data.coeff_seper, 1) ~= Nmu
  error('validate_coeff_coarse_vs_fine:shape', ...
    'coeff_seper row count (%d) mismatches R_sampling_RLU rows (%d).', ...
    size(isdf_data.coeff_seper, 1), Nmu);
end

if size(wf_data.c, 2) ~= nb
  error('validate_coeff_coarse_vs_fine:band', ...
    'wf_data.c nb (%d) mismatches coeff_seper nb (%d).', ...
    size(wf_data.c, 2), nb);
end

tol_integer = 1e-10;
R_round = round(R_sampling_RLU);
if max(abs(R_sampling_RLU(:) - R_round(:))) > tol_integer
  error('validate_coeff_coarse_vs_fine:nonIntegerSampling', ...
    ['R_sampling_RLU has non-integer entries, cannot map to exact fine-grid indices. ' ...
     'max deviation = %.3e'], ...
    max(abs(R_sampling_RLU(:) - R_round(:))));
end
R_int = int32(R_round);

nx = fftgrid(1);
ny = fftgrid(2);
nz = fftgrid(3);

abs_err = zeros(Nmu, 1);
rel_err = zeros(Nmu, 1);
ir_map = zeros(Nmu, 1, 'int32');

for imu = 1:Nmu
  r = mod(R_int(imu, :) + fftgrid, fftgrid);
  ir = int32(1 + r(1) + r(2) * nx + r(3) * nx * ny);
  if ir < 1 || ir > nx * ny * nz
    error('validate_coeff_coarse_vs_fine:irRange', ...
      'Mapped ir out of range at imu=%d.', imu);
  end
  ir_map(imu) = ir;

  c_coarse = squeeze(isdf_data.coeff_seper(imu, :, 1, 1));
  c_fine = squeeze(wf_data.c(ir, :, 1, 1));
  diff_vec = c_coarse(:) - c_fine(:);
  abs_err(imu) = norm(diff_vec);
  rel_err(imu) = abs_err(imu) / max(norm(c_fine(:)), eps);
end

[max_abs_err, imu_max_abs] = max(abs_err);
[max_rel_err, imu_max_rel] = max(rel_err);

fprintf('\n=== validate_coeff_coarse_vs_fine ===\n');
fprintf('id_nn = %d\n', id_nn);
fprintf('fine fftgrid   = [%d %d %d]\n', fftgrid(1), fftgrid(2), fftgrid(3));
fprintf('coarse fftgrid = [%d %d %d]\n', fftgrid_c(1), fftgrid_c(2), fftgrid_c(3));
fprintf('Nmu = %d, nb = %d\n', Nmu, nb);
fprintf('max |c_coarse-c_fine|_2 = %.6e (imu=%d, ir=%d)\n', ...
  max_abs_err, imu_max_abs, ir_map(imu_max_abs));
fprintf('max rel err           = %.6e (imu=%d, ir=%d)\n', ...
  max_rel_err, imu_max_rel, ir_map(imu_max_rel));
fprintf('mean abs err          = %.6e\n', mean(abs_err));
fprintf('mean rel err          = %.6e\n', mean(rel_err));
fprintf('=== end ===\n\n');

