% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20

function root = adaptive_checkpoint_storage_root()
%ADAPTIVE_CHECKPOINT_STORAGE_ROOT  Directory for isdf_adaptive_checkpoint_*.mat files.

  cfile = fullfile(pwd, 'SAVE', 'config.mat');
  if exist(cfile, 'file') == 2
    S = load(cfile, 'config');
    if isfield(S, 'config') && isfield(S.config, 'CONTROL') && ...
        isfield(S.config.CONTROL, 'storage_dir')
      r = strtrim(char(string(S.config.CONTROL.storage_dir)));
      if ~isempty(r)
        root = r;
        return
      end
    end
  end
  root = fullfile(pwd, 'SAVE');
end
