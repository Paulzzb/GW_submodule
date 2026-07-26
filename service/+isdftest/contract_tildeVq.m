function val = contract_tildeVq(id, iqibz, c)
% CONTRACT_TILDEVQ  val = c' tildeVq c with c already transformed (from get_rho_xalpha).
%
%   val = isdftest.contract_tildeVq(id, iqibz, c)
%
% Prefer passing c from isdftest.get_rho_xalpha (already c_t = C^{-1/2} c).

  isdf_data = isdftest.get(id);
  if isempty(isdf_data.tildeVq)
    error('isdftest:contract_tildeVq:MissingData', ...
      'tildeVq empty for id=%d; run isdftest.gen_tildeVq first.', int32(id));
  end
  c = c(:);
  tildeVq = isdf_data.tildeVq(:, :, iqibz);
  if size(tildeVq, 1) ~= numel(c)
    error('isdftest:contract_tildeVq:Size', ...
      'c length %d inconsistent with tildeVq (%d).', numel(c), size(tildeVq, 1));
  end
  val = c' * tildeVq * c;
end
