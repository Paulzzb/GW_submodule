function tildeWq = gen_tildeWq(id_vc, iqibz, Kq, id_outer, flagherm)
% Kq is calculated by gen_Kq.m, which is using vc data
% Now, formally,
%     tildeW_q = <outer|Vq|in> * Kq^{-1} * <in|Vq|outer> + <outer|Vq|outer>,
% where <outer|Vq|outer> is calculated by gen_tildeVq in isdf.get(id_outer),
% and <outer|Vq|in> will be calculated here.
% 
% Normally, we only need to calculate id_outer.desc = 'vn'. However, in case
% of COHSEX approximation, we also need to calculate id_outer.desc = 'nn'.

if nargin < 4
  error('gen_tildeWq: Missing input: id_vc, iqibz, Kq, id_outer');
end

if nargin < 5
  flagherm = false;
end

vc_data = isdf.get(id_vc);
outer_data = isdf.get(id_outer);
%
vc_fac = vc_data.CCHq_trunc_factors{double(iqibz)};
Nkeep_vc = vc_fac.N_keep;
outer_fac = outer_data.CCHq_trunc_factors{double(iqibz)};
Nkeep_outer = outer_fac.N_keep;
%
coulomb_data = coulomb.get();
%
vcoul_q = coulomb_data.vcoul(:, iqibz);
if iqibz == 1
  vcoul_q(1) = coulomb_data.vcoul0;
end

% vcoul_q = vcoul_q * 13.6059;
%
tildeWq = zeros(Nkeep_outer, Nkeep_outer);
% Calculate <outer|Vq|in>
helperqG_vc = vc_data.helperqG(:, :, iqibz);
helperqG_outer = outer_data.helperqG(:, :, iqibz);

vc_Vq_outer = helperqG_vc(:, 1:Nkeep_vc)' * diag(vcoul_q) * helperqG_outer(:, 1:Nkeep_outer);

%
% t = vc_Vq_outer' * inv(Kq) * vc_Vq_outer;

if flagherm
  try
    L_Kq = chol(Kq, "lower");
    vc_tmp = L_Kq \ vc_Vq_outer;
    tildeWq = vc_tmp' * vc_tmp;
  catch
    tildeWq = vc_Vq_outer' * (Kq \ vc_Vq_outer);
  end
else
  tildeWq = vc_Vq_outer' * (Kq \ vc_Vq_outer);
end

% if norm(t - tildeWq, 'fro') / norm(t, 'fro') > 1e-6
%   warning('gen_tildeWq: inconsistent t and tildeWq, rel=%.6e', norm(t - tildeWq, 'fro') / norm(t, 'fro'));
% end
% tildeWq = tildeWq + outer_data.tildeVq(:, :, iqibz);



end