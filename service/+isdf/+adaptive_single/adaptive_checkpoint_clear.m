% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20

function adaptive_checkpoint_clear(coarse_id, desc)
%ADAPTIVE_CHECKPOINT_CLEAR  Remove checkpoint file(s) for a coarse ISDF slot.
%
%   isdf.adaptive_single.adaptive_checkpoint_clear(coarse_id)
%   isdf.adaptive_single.adaptive_checkpoint_clear(coarse_id, desc)
%
% Deletes desc-keyed checkpoint and legacy id-only file (if present).

  if nargin < 2 || isempty(desc)
    coarse_data = isdf.get(coarse_id);
    desc = char(coarse_data.desc);
  else
    desc = char(string(desc));
  end

  paths = { ...
    isdf.adaptive_single.adaptive_checkpoint_path(coarse_id, desc), ...
    isdf.adaptive_single.adaptive_checkpoint_legacy_path(coarse_id) ...
  };
  for k = 1:numel(paths)
    if exist(paths{k}, 'file') == 2
      delete(paths{k});
      fprintf(1, 'adaptiveisdf: removed checkpoint %s\n', paths{k});
    end
  end
end
