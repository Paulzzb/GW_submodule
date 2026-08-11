% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/11 ZZ

function clean_case_outputs(case_dir)
%CLEAN_CASE_OUTPUTS  Remove prior computed outputs; keep namelist + groundstate.
%
%   clean_case_outputs(CASE_DIR)
%
%   Deletes CASE_DIR/SAVE/, CASE_DIR/isdf_report/, and root artifacts
%   (logs / qp*). Leaves ./test intact.
%
%   Safe before QPstartup: adds util/ only if filename_map is not on the path.

  examples_dir = fileparts(mfilename('fullpath'));
  gw_root = fileparts(examples_dir);
  if isempty(which('filename_map'))
    addpath(fullfile(gw_root, 'util'));
  end

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

  patterns = { ...
    'qp.dat', ...
    'qp_*.dat', ...
    'r-*.log', ...
    'l-*.log', ...
    };
  for p = 1:numel(patterns)
    hits = dir(fullfile(case_dir, patterns{p}));
    for k = 1:numel(hits)
      name = hits(k).name;
      if strcmp(name, '.') || strcmp(name, '..') || hits(k).isdir
        continue
      end
      if isfield(hits, 'folder') && ~isempty(hits(k).folder)
        f = fullfile(hits(k).folder, name);
      else
        f = fullfile(case_dir, name);
      end
      if ~isfile(f)
        continue
      end
      fprintf('clean: removing %s\n', f);
      delete(f);
    end
  end
end
