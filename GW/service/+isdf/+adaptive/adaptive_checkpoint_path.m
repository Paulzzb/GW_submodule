% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/25

function fpath = adaptive_checkpoint_path(coarse_id)
%ADAPTIVE_CHECKPOINT_PATH  Full path to adaptive ISDF checkpoint for coarse pool id.
%
%   fpath = isdf.adaptive.adaptive_checkpoint_path(coarse_id)
%
% File name: isdf_adaptive_checkpoint_id<coarse_id>.mat under storage root.
% Storage root: config.CONTROL.storage_dir from ./SAVE/config.mat if present,
% otherwise <pwd>/SAVE.

  root = adaptive_checkpoint_storage_root();
  fpath = fullfile(root, sprintf('isdf_adaptive_checkpoint_id%d.mat', int32(coarse_id)));
end

function root = adaptive_checkpoint_storage_root()
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
