% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20

function adaptive_checkpoint_save(coarse_id, id_adaptive)
%ADAPTIVE_CHECKPOINT_SAVE  Write adaptive ISDF pool slot to checkpoint .mat.
%
%   isdf.adaptive_single.adaptive_checkpoint_save(coarse_id, id_adaptive)
%
% Saves isdf_data, id_coarse, id_adaptive, and isdf_desc (vc/vn/nn of coarse slot).

  coarse_data = isdf.get(coarse_id);
  isdf_desc = char(coarse_data.desc);
  fpath = isdf.adaptive_single.adaptive_checkpoint_path(coarse_id, isdf_desc);
  isdf_data = isdf.get(id_adaptive);
  id_coarse = coarse_id;
  parentDir = fileparts(fpath);
  if ~isempty(parentDir) && exist(parentDir, 'dir') ~= 7
    mkdir(parentDir);
  end
  save(fpath, 'isdf_data', 'id_coarse', 'id_adaptive', 'isdf_desc', '-v7.3');
  fprintf(1, ['adaptiveisdf: saved checkpoint %s (id_coarse=%d, id_adaptive=%d, ', ...
    'desc=%s).\n'], fpath, double(coarse_id), double(id_adaptive), isdf_desc);
end
