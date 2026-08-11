% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20

function fpath = adaptive_checkpoint_path(coarse_id, desc)
%ADAPTIVE_CHECKPOINT_PATH  Checkpoint path under storage_dir.
%
%   fpath = isdf.adaptive.adaptive_checkpoint_path(coarse_id)
%   fpath = isdf.adaptive.adaptive_checkpoint_path(coarse_id, desc)
%   isdf.adaptive.adaptive_checkpoint_path('set', storage_dir)
%   isdf.adaptive.adaptive_checkpoint_path('clear')
%
% File: <storage_dir>/isdf_adaptive_checkpoint_<desc>_id<coarse_id>.mat

  persistent storage_dir

  if nargin >= 1 && (ischar(coarse_id) || isstring(coarse_id))
    switch lower(char(string(coarse_id)))
      case 'set'
        storage_dir = char(string(desc));
        fpath = storage_dir;
        return
      case 'clear'
        storage_dir = [];
        fpath = '';
        return
    end
  end

  if nargin < 2 || isempty(desc)
    coarse_data = isdf.get(coarse_id);
    desc = char(coarse_data.desc);
  else
    desc = char(string(desc));
  end

  if ~isempty(storage_dir)
    root = storage_dir;
  else
    root = default_param_values().CONTROL.storage_dir;
  end
  fpath = fullfile(root, sprintf('isdf_adaptive_checkpoint_%s_id%d.mat', desc, int32(coarse_id)));
end
