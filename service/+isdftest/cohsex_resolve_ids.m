function [id_vc, id_vn, id_nn] = cohsex_resolve_ids(config)
%COHSEX_RESOLVE_IDS  Map GW config to ISDF pool slots for static COHSEX (service +isdf).
%
%   [id_vc, id_nn] = isdf.cohsex_resolve_ids(config)
%
% Resolution order for each label ('vc', 'vn', 'nn'):
%   1) config.ISDF.id_vc / id_vn / id_nn if nonempty
%   2) First assigned pool entry whose desc matches (case-sensitive char)
%
% Requires an initialized isdftest.manager pool (relay / driver already ran).

  id_vc = local_pick(config, 'id_vc', "vc");
  id_vn = local_pick(config, 'id_vn', "vn");
  id_nn = local_pick(config, 'id_nn', "nn");
  if isempty(id_vc) && isempty(id_vn) && isempty(id_nn)
    error('isdf:cohsex_resolve_ids', 'No ISDF slot assigned for any label.');
  end
end

function id = local_pick(~, ~, want_desc)
  L = isdftest.manager('list');
  for k = 1:numel(L)
    if L(k).empty || ~L(k).assigned
      continue
    end
    if strcmp(char(L(k).desc), want_desc)
      id = int32(L(k).id);
      return
    end
  end
  id = [];
end
