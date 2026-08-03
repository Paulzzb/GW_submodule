% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20

function fpath = adaptive_checkpoint_legacy_path(coarse_id)
%ADAPTIVE_CHECKPOINT_LEGACY_PATH  Old checkpoint name (id only, no desc).

  root = isdf.adaptive.adaptive_checkpoint_storage_root();
  fpath = fullfile(root, sprintf('isdf_adaptive_checkpoint_id%d.mat', int32(coarse_id)));
end
