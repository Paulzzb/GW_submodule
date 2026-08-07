% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/07/31 ZZ

function clean_case_outputs(case_dir)
%CLEAN_CASE_OUTPUTS  Remove prior computed outputs; keep namelist + groundstate.
%
%   clean_case_outputs(CASE_DIR)
%
%   Deletes CASE_DIR/SAVE/ and common run artifacts next to the namelist.
%   Leaves ./test and ./qe.save intact.
%
%   Safe before QPstartup: adds util_profile so filename_map is visible.

  tests_dir = fileparts(mfilename('fullpath'));
  gw_root = fileparts(tests_dir);
  addpath(fullfile(gw_root, 'util_profile'));

  save_dir = fullfile(case_dir, 'SAVE');
  if isfolder(save_dir)
    fprintf('clean: removing %s\n', save_dir);
    rmdir(save_dir, 's');
  end

  def = filename_map();
  report_dir = fullfile(case_dir, def.isdf_report_dir);
  if isfolder(report_dir)
    fprintf('clean: removing %s\n', report_dir);
    rmdir(report_dir, 's');
  end

  % Legacy root-level OF names (pre-isdf_report_dir layout) + shared artifacts.
  hf_glob = strrep(def.hf_report, '%d', '*');
  ad_glob = strrep(def.adaptive_report, '%d', '*');
  patterns = { ...
    'adaptiveisdf_*.txt', ...
    'isdf_validate_*.txt', ...
    hf_glob, ...
    ad_glob, ...
    def.cond_report, ...
    [def.cond_report '_*'], ...
    'qp_*.dat', ...
    '*_relay_stage.mat', ...
    'r-*.log', ...
    'l-*.log', ...
    'Ex.dat' ...
    };
  for p = 1:numel(patterns)
    hits = dir(fullfile(case_dir, patterns{p}));
    for k = 1:numel(hits)
      if hits(k).isdir
        continue
      end
      f = fullfile(case_dir, hits(k).name);
      fprintf('clean: removing %s\n', f);
      delete(f);
    end
  end
end
