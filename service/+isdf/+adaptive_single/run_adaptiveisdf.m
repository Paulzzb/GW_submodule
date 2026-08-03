% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/16

function run_adaptiveisdf(type, output_dir)
% run_adaptiveisdf  Run adaptive ISDF from a relay profile (same idea as test_*/run_adaptiveisdf_smoke).
%
%   isdf.adaptive_single.run_adaptiveisdf('nn', output_dir)
%   isdf.adaptive_single.run_adaptiveisdf('vn', output_dir)
%
% Inputs:
%   type       鈥?'nn' or 'vn' (case-insensitive). Selects coarse ISDF slot description and nisdf rule.
%   output_dir 鈥?Directory where isdf_validate_HF writes the HF report via isdf.report.gen_report
%                (isdf_validate_HF_id<idnew>.txt for the new adaptive pool id).
%
% Usage:
%   Set MATLAB current folder to a GW test profile that contains SAVE/config.mat and
%   test_relay_stage.mat (e.g. test_Si), then call isdf.adaptive_single.run_adaptiveisdf(...).
%
% This function temporarily cd's to output_dir only while adaptiveisdf runs so that
% validation uses that directory as pwd for report generation.

  profile_dir = pwd;

  type = lower(strtrim(char(string(type))));
  if ~ismember(type, {'vn', 'nn'})
    error('run_adaptiveisdf:type', 'type must be ''vn'' or ''nn'' (got ''%s'').', type);
  end

  output_dir = char(string(output_dir));
  if isempty(strtrim(output_dir))
    error('run_adaptiveisdf:output_dir', 'output_dir must be a non-empty path.');
  end
  if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
  end

  cfg_path = fullfile(profile_dir, 'SAVE', 'config.mat');
  if ~exist(cfg_path, 'file')
    error('run_adaptiveisdf:config', 'Missing %s (expected under current folder %s).', ...
      cfg_path, profile_dir);
  end
  load(cfg_path, 'config');

  if ~isfield(config, 'ISDF') || ~config.ISDF.isisdf
    error('run_adaptiveisdf:isisdf', 'Set isisdf in test input for this run.');
  end

  stage_path = fullfile(profile_dir, 'test_relay_stage.mat');
  if ~exist(stage_path, 'file')
    error('run_adaptiveisdf:stage', 'Missing %s. Run input_driver first.', stage_path);
  end

  service_reset_persistent();
  relay.stage_from_db(stage_path);
  relay.restore();

  isdf.debug.init_from_config(config);

  cfg = config.ISDF;

  isdf.free();

  if strcmp(type, 'vn')
    desc_token = 'vn';
    id_slot = isdf.isdf_add('vn');
    isdf.set_nrange(id_slot, config.SYSTEM);
    slot_data = isdf.get(id_slot);
    nmu_target = cfg.isdf_ratio_type2 * sqrt(double(length(slot_data.nrange1)) * double(length(slot_data.nrange2)));
    slot_data.nisdf = int32(max(1, ceil(nmu_target)));
    slot_data.assigned = true;
    isdf.save2mod(slot_data, id_slot);
  else
    desc_token = 'nn';
    id_slot = isdf.isdf_add('nn');
    isdf.set_nrange(id_slot, config.SYSTEM);
    slot_data = isdf.get(id_slot);
    nmu_target = cfg.isdf_ratio_type3 * sqrt(double(length(slot_data.nrange1)) * double(length(slot_data.nrange2)));
    slot_data.nisdf = int32(max(1, ceil(nmu_target)));
    slot_data.assigned = true;
    isdf.save2mod(slot_data, id_slot);
  end

  isdf.coeff.gen_coeff(cfg, id_slot);
  isdf.coeff.print_coarse_grid_report(id_slot);

  L = isdf.manager('list');
  id_slot = [];
  for k = 1:numel(L)
    if L(k).assigned && strcmp(char(L(k).desc), desc_token)
      id_slot = L(k).id;
      break;
    end
  end
  if isempty(id_slot)
    error('run_adaptiveisdf:id', 'No assigned ISDF slot with desc ''%s''.', desc_token);
  end

  isdf.adaptive_single.isdf_schur_update('clear');
  isdf.adaptive_single.adaptive_weight('clear');

  cleanup_cd = onCleanup(@() cd(profile_dir));
  cd(output_dir);
  isdf.adaptive_single.adaptiveisdf(id_slot, cfg);

  rep_hf = dir(fullfile(output_dir, 'isdf_validate_HF_id*.txt'));
  if isempty(rep_hf)
    warning('run_adaptiveisdf:noHFReport', ...
      ['No isdf_validate_HF_id*.txt in output_dir=%s. ', ...
       'Check adaptiveisdf / isdf_validation errors.'], output_dir);
  else
    fprintf('run_adaptiveisdf: HF report(s) in output_dir:\n');
    for ri = 1:numel(rep_hf)
      fprintf('  %s\n', fullfile(output_dir, rep_hf(ri).name));
    end
  end

  rep_ad = dir(fullfile(output_dir, sprintf('adaptiveisdf_id%d.txt', id_slot)));
  if isempty(rep_ad)
    warning('run_adaptiveisdf:noAdaptiveReport', ...
      'Expected adaptiveisdf_id%d.txt in output_dir=%s.', id_slot, output_dir);
  else
    fprintf('run_adaptiveisdf: adaptiveisdf phase-1 report:\n');
    fprintf('  %s\n', fullfile(output_dir, rep_ad(1).name));
  end

  fprintf('run_adaptiveisdf: OK (type=%s, coarse id=%d, profile_dir=%s)\n', ...
    type, id_slot, profile_dir);
end
