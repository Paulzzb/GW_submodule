% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/10 ZZ

function clean_case_outputs(case_dir)
%CLEAN_CASE_OUTPUTS  Remove prior computed outputs; keep namelist + groundstate.
%
%   clean_case_outputs(CASE_DIR)
%
%   Deletes CASE_DIR/SAVE/, CASE_DIR/isdf_report/, and shared root artifacts
%   (logs / qp / Ex / relay). Leaves ./test and ./qe.save intact.
%
%   If CASE_DIR/SAVE is a symlink (shared hub SAVE), only the link is left
%   alone — never rmdir into the target (would wipe Si_gamma/SAVE).
%
%   Safe before QPstartup: adds util/ so filename_map is visible.

  tests_dir = fileparts(mfilename('fullpath'));
  gw_root = fileparts(tests_dir);
  addpath(fullfile(gw_root, 'util'));

  save_dir = fullfile(case_dir, 'SAVE');
  if local_is_symlink(save_dir)
    fprintf('clean: keep SAVE symlink %s\n', save_dir);
  elseif isfolder(save_dir)
    fprintf('clean: removing %s\n', save_dir);
    rmdir(save_dir, 's');
  end

  def = filename_map();
  report_dir = fullfile(case_dir, def.isdf_report_dir);
  if isfolder(report_dir)
    fprintf('clean: removing %s\n', report_dir);
    rmdir(report_dir, 's');
  end

  % Root-level artifacts only (ISDF OF live under isdf_report/, already removed).
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

function tf = local_is_symlink(p)
  tf = false;
  try
    tf = java.nio.file.Files.isSymbolicLink(java.nio.file.Paths.get(p));
  catch
    if isunix
      [st, ~] = system(sprintf('test -L %s', local_shell_quote(p)));
      tf = (st == 0);
    end
  end
end

function s = local_shell_quote(p)
  s = ['''' strrep(p, '''', '''\'''''') ''''];
end
