% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/02

function [loaded, idnew] = adaptive_checkpoint_import_from_isdf(coarse_id, source)
%ADAPTIVE_CHECKPOINT_IMPORT_FROM_ISDF  Load +isdf adaptive checkpoint into +isdftest pool.
%
%   [loaded, idnew] = isdftest.adaptive.adaptive_checkpoint_import_from_isdf(coarse_id, source)
%
% source: directory containing isdf_adaptive_checkpoint_<desc>_id<N>.mat, or a full .mat path.
% On success, writes a native +isdftest checkpoint under the current storage root.

  loaded = false;
  idnew = [];

  coarse_data = isdftest.get(coarse_id);
  expected_desc = char(coarse_data.desc);

  fpath = local_resolve_isdf_checkpoint_path(source, coarse_id, expected_desc);
  if isempty(fpath) || exist(fpath, 'file') ~= 2
    return;
  end

  S = load(fpath);
  if ~isfield(S, 'isdf_data')
    error('isdftest:adaptive:importIsdf:BadFile', ...
      'Checkpoint %s must contain variable isdf_data.', fpath);
  end
  if ~isa(S.isdf_data, 'isdf.base.isdf_m')
    error('isdftest:adaptive:importIsdf:Type', ...
      'Checkpoint isdf_data must be isdf.base.isdf_m (got %s).', class(S.isdf_data));
  end

  ck_desc = '';
  if isfield(S, 'isdf_desc')
    ck_desc = char(string(S.isdf_desc));
  elseif isprop(S.isdf_data, 'desc')
    ck_desc = char(S.isdf_data.desc);
  end
  if ~isempty(ck_desc) && ~strcmp(ck_desc, expected_desc)
    fprintf(1, ['adaptiveisdf: skip +isdf import %s (isdf_desc=''%s'', expected ''%s'').\n'], ...
      fpath, ck_desc, expected_desc);
    return;
  end

  d = isdftest.adaptive.convert_isdf_m_to_isdftest_m(S.isdf_data);
  idnew = isdftest.isdftest_add('adaptive');
  d.id = idnew;
  isdftest.save2mod(d, idnew);
  loaded = true;

  isdftest.adaptive.adaptive_checkpoint_save(coarse_id, idnew);

  if isfield(S, 'id_coarse') && isfield(S, 'id_adaptive')
    fprintf(1, ['adaptiveisdf: imported +isdf checkpoint %s -> isdftest pool id=%d ', ...
      '(file id_coarse=%s, id_adaptive=%s, desc=%s).\n'], ...
      fpath, double(idnew), mat2str(S.id_coarse), mat2str(S.id_adaptive), expected_desc);
  else
    fprintf(1, 'adaptiveisdf: imported +isdf checkpoint %s -> isdftest pool id=%d (desc=%s).\n', ...
      fpath, double(idnew), expected_desc);
  end
end

function fpath = local_resolve_isdf_checkpoint_path(source, coarse_id, expected_desc)
  fpath = '';
  if nargin < 1 || isempty(source)
    return;
  end
  source = char(string(source));
  if exist(source, 'file') == 2
    fpath = source;
    return;
  end
  if exist(source, 'dir') ~= 7
    return;
  end
  fpath = fullfile(source, sprintf('isdf_adaptive_checkpoint_%s_id%d.mat', ...
    expected_desc, int32(coarse_id)));
  if exist(fpath, 'file') == 2
    return;
  end
  legacy = fullfile(source, sprintf('isdf_adaptive_checkpoint_id%d.mat', int32(coarse_id)));
  if exist(legacy, 'file') == 2 && local_legacy_desc_matches(legacy, expected_desc)
    fpath = legacy;
  else
    fpath = '';
  end
end

function ok = local_legacy_desc_matches(legacy_path, expected_desc)
  ok = false;
  info = whos('-file', legacy_path);
  names = {info.name};
  if ~any(strcmp(names, 'isdf_desc'))
    return;
  end
  S = load(legacy_path, 'isdf_desc');
  ok = strcmp(char(string(S.isdf_desc)), expected_desc);
end
