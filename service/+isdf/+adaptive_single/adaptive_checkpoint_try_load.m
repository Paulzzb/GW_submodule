% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20

function [loaded, idnew] = adaptive_checkpoint_try_load(coarse_id)
%ADAPTIVE_CHECKPOINT_TRY_LOAD  Restore adaptive ISDF from checkpoint if file exists.
%
%   [loaded, idnew] = isdf.adaptive_single.adaptive_checkpoint_try_load(coarse_id)
%
% Checkpoint must match the coarse slot desc (vc/vn/nn). Legacy id-only files
% are used only when isdf_desc in the file matches the current coarse slot.

  loaded = false;
  idnew = [];

  coarse_data = isdf.get(coarse_id);
  expected_desc = char(coarse_data.desc);

  fpath = isdf.adaptive_single.adaptive_checkpoint_path(coarse_id, expected_desc);
  if exist(fpath, 'file') ~= 2
    legacy = isdf.adaptive_single.adaptive_checkpoint_legacy_path(coarse_id);
    if exist(legacy, 'file') ~= 2
      return
    end
    if ~local_legacy_desc_matches(legacy, expected_desc)
      fprintf(1, ['adaptiveisdf: skip legacy checkpoint %s (desc mismatch or ', ...
        'missing isdf_desc; expected ''%s'').\n'], legacy, expected_desc);
      return
    end
    fpath = legacy;
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

  ck_desc = '';
  if isfield(S, 'isdf_desc')
    ck_desc = char(string(S.isdf_desc));
  elseif isprop(S.isdf_data, 'desc')
    ck_desc = char(S.isdf_data.desc);
  end
  if ~isempty(ck_desc) && ~strcmp(ck_desc, expected_desc)
    fprintf(1, ['adaptiveisdf: skip checkpoint %s (isdf_desc=''%s'', expected ''%s'').\n'], ...
      fpath, ck_desc, expected_desc);
    return
  end

  idnew = isdf.isdf_add('adaptive');
  d = S.isdf_data;
  d.id = idnew;
  isdf.save2mod(d, idnew);
  loaded = true;

  if isfield(S, 'id_coarse') && isfield(S, 'id_adaptive')
    fprintf(1, ['adaptiveisdf: loaded checkpoint %s -> pool id=%d ', ...
      '(file id_coarse=%s, id_adaptive=%s, desc=%s).\n'], ...
      fpath, double(idnew), mat2str(S.id_coarse), mat2str(S.id_adaptive), expected_desc);
  else
    fprintf(1, 'adaptiveisdf: loaded checkpoint %s -> pool id=%d (desc=%s).\n', ...
      fpath, double(idnew), expected_desc);
  end
end

function ok = local_legacy_desc_matches(legacy_path, expected_desc)
  ok = false;
  info = whos('-file', legacy_path);
  names = {info.name};
  if any(strcmp(names, 'isdf_desc'))
    S = load(legacy_path, 'isdf_desc');
    ok = strcmp(char(string(S.isdf_desc)), expected_desc);
    return
  end
  % Legacy files without isdf_desc: do not load (may be wrong type on same pool id).
end
