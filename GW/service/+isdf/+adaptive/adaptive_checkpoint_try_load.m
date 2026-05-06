% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/25

function [loaded, idnew] = adaptive_checkpoint_try_load(coarse_id)
%ADAPTIVE_CHECKPOINT_TRY_LOAD  Restore adaptive ISDF from checkpoint if file exists.
%
%   [loaded, idnew] = isdf.adaptive.adaptive_checkpoint_try_load(coarse_id)
%
% If checkpoint file exists: allocates new 'adaptive' pool slot, save2mod loaded
% isdf_data, returns loaded=true and idnew. Otherwise loaded=false, idnew=[].

  loaded = false;
  idnew = [];
  fpath = isdf.adaptive.adaptive_checkpoint_path(coarse_id);
  if exist(fpath, 'file') ~= 2
    return
  end

  S = load(fpath);
  if ~isfield(S, 'isdf_data')
    error('isdf:adaptive:checkpoint:BadFile', ...
      'Checkpoint %s must contain variable isdf_data (isdf.base.isdf_m).', fpath);
  end
  if ~isa(S.isdf_data, 'isdf.base.isdf_m')
    error('isdf:adaptive:checkpoint:Type', ...
      'Checkpoint isdf_data must be isdf.base.isdf_m (got %s).', class(S.isdf_data));
  end

  idnew = isdf.isdf_add('adaptive');
  d = S.isdf_data;
  d.id = idnew;
  isdf.save2mod(d, idnew);
  loaded = true;

  if isfield(S, 'id_coarse') && isfield(S, 'id_adaptive')
    fprintf(1, ['adaptiveisdf: loaded checkpoint %s -> pool id=%d ', ...
      '(file id_coarse=%s, id_adaptive=%s).\n'], ...
      fpath, double(idnew), mat2str(S.id_coarse), mat2str(S.id_adaptive));
  else
    fprintf(1, 'adaptiveisdf: loaded checkpoint %s -> pool id=%d.\n', fpath, double(idnew));
  end
end
