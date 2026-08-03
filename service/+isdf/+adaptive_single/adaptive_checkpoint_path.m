% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20

function fpath = adaptive_checkpoint_path(coarse_id, desc)
%ADAPTIVE_CHECKPOINT_PATH  Checkpoint path keyed by ISDF desc + coarse pool id.
%
%   fpath = isdf.adaptive_single.adaptive_checkpoint_path(coarse_id)
%   fpath = isdf.adaptive_single.adaptive_checkpoint_path(coarse_id, desc)
%
% File name: isdf_adaptive_checkpoint_<desc>_id<coarse_id>.mat under storage root.
% (Legacy name isdf_adaptive_checkpoint_id<coarse_id>.mat is read-only fallback in try_load.)

  if nargin < 2 || isempty(desc)
    coarse_data = isdf.get(coarse_id);
    desc = char(coarse_data.desc);
  else
    desc = char(string(desc));
  end

  root = isdf.adaptive_single.adaptive_checkpoint_storage_root();
  fpath = fullfile(root, sprintf('isdf_adaptive_checkpoint_%s_id%d.mat', desc, int32(coarse_id)));
end
