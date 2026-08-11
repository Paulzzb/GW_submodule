% Demo entry:
% 1) Restore service persistent data from staged snapshot.
% 2) Build k1-k2 orbit representatives with shared G0 and unordered pairs.
% 3) Save artifacts for inspection.

this_dir = '';
% Prefer current folder when running from test_GaN directly.
if exist(fullfile(pwd, 'SAVE'), 'dir') && exist(fullfile(pwd, 'test'), 'file')
  this_dir = pwd;
end
if isempty(this_dir)
  this_file = mfilename('fullpath');
  if isempty(this_file)
    this_file = which(mfilename);
  end
  this_dir = fileparts(this_file);
end
cleanup_return_dir = onCleanup(@() cd(this_dir)); %#ok<NASGU>
service_dir = fullfile(this_dir, '..', '..', 'service');
pair_sym_dir = fullfile(service_dir, '+pair_symmetry');
if ~exist(pair_sym_dir, 'dir')
  error('demo:MissingPairSymmetryDir', 'Cannot find +pair_symmetry directory: %s', pair_sym_dir);
end
addpath(service_dir);
cleanup_obj = onCleanup(@() rmpath(service_dir)); %#ok<NASGU>
rehash;

if exist('pair_symmetry.driver', 'file') ~= 2
  orig_dir = pwd;
  cleanup_cd = onCleanup(@() cd(orig_dir)); %#ok<NASGU>
  cd(service_dir);
end

cd(fullfile(this_dir, '..', '..'));
QPstartup;
cd(this_dir);

service_reset_persistent;
packages_reset_persistent;
test_stage;

pair_symmetry.driver([], struct());
pair_data = pair_symmetry.get();

k1k2_representation = pair_data.representation;
k1k2_mapping = pair_data.mapping;
marked = pair_data.mapping(:, :, 2) > 0;
tol = 1e-5;
validation_report = struct();
validation_report.ok = true;
validation_report.note = 'Validation is skipped in package-mode demo.';

save_path = fullfile(this_dir, 'SAVE', 'k1k2_representation.mat');


pair_symmetry.save_k1k2_result(save_path, ...
  k1k2_representation, k1k2_mapping, marked, tol, pair_data.nk, int32(0), validation_report);
