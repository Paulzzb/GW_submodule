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

  save_dir = fullfile(case_dir, 'SAVE');
  if isfolder(save_dir)
    fprintf('clean: removing %s\n', save_dir);
    rmdir(save_dir, 's');
  end

  patterns = { ...
    'adaptiveisdf_*.txt', ...
    'isdf_validate_*.txt', ...
    'o-ISDF_HF_id*', ...
    'o-ISDF_adaptive_id*', ...
    'o-ISDF_cond', ...
    'o-ISDF_cond_*', ...
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
