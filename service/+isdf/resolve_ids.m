% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06 ZZ

function [id_vc, id_vn, id_nn] = resolve_ids(config)
%RESOLVE_IDS  Map GW config to ISDF pool slots (vc / vn / nn).
%
%   [id_vc, id_vn, id_nn] = isdf.resolve_ids(config)
%
% Resolution order for each label ('vc', 'vn', 'nn'):
%   1) config.ISDF.id_vc / id_vn / id_nn if nonempty
%   2) First assigned pool entry whose desc matches (case-sensitive char)
%
% Requires an initialized isdf.manager pool (relay / driver already ran).

  id_vc = local_pick(config, 'id_vc', "vc");
  id_vn = local_pick(config, 'id_vn', "vn");
  id_nn = local_pick(config, 'id_nn', "nn");
  if isempty(id_vc) && isempty(id_vn) && isempty(id_nn)
    error('isdf:resolve_ids', 'No ISDF slot assigned for any label.');
  end
end

function id = local_pick(~, ~, want_desc)
  L = isdf.manager('list');
  for k = 1:numel(L)
    if L(k).empty || ~L(k).assigned
      continue
    end
    d = char(string(L(k).desc));
    if strcmp(d, want_desc) || startsWith(d, [char(want_desc) '_'])
      id = int32(L(k).id);
      return
    end
  end
  id = [];
end
