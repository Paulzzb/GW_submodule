function val = contract_tildeVq(id, iqibz, c)
% CONTRACT_TILDEVQ  val = c' tildeVq c with c already transformed (from get_rho_xalpha).
%
%   val = isdf.contract_tildeVq(id, iqibz, c)
%
% Prefer passing c from isdf.get_rho_xalpha (already c_t = C^{-1/2} c).

  isdf_data = isdf.get(id);
  if isempty(isdf_data.tildeVq)
    error('isdf:contract_tildeVq:MissingData', ...
      'tildeVq empty for id=%d; run isdf.gen_tildeVq first.', int32(id));
  end
  c = c(:);
  tildeVq = isdf_data.tildeVq(:, :, iqibz);
  if size(tildeVq, 1) ~= numel(c)
    error('isdf:contract_tildeVq:Size', ...
      'c length %d inconsistent with tildeVq (%d).', numel(c), size(tildeVq, 1));
  end
  val = c' * tildeVq * c;
end
