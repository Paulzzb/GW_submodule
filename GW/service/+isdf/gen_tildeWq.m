function tildeWq = gen_tildeWq(id_vc, iqibz, Kq, id_outer)
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

vc_data = isdf.get(id_vc);
outer_data = isdf.get(id_outer);
%
Nisdf_vc = vc_data.nisdf;
Nisdf_outer = outer_data.nisdf;
%
coulomb_data = coulomb.get();
%
vcoul_q = coulomb_data.vcoul(:, iqibz);
if iqibz == 1
  vcoul_q(1) = coulomb_data.vcoul0;
end

% vcoul_q = vcoul_q * 13.6059;
%
tildeWq = zeros(Nisdf_outer, Nisdf_outer);
% Calculate <outer|Vq|in>
helperqG_vc = vc_data.helperqG(:, :, iqibz);
helperqG_outer = outer_data.helperqG(:, :, iqibz);
vc_Vq_outer = helperqG_vc' * diag(vcoul_q) * helperqG_outer;
%
% t = vc_Vq_outer' * inv(Kq) * vc_Vq_outer;
L_Kq = chol(Kq, "lower");
vc_Vq_outer = L_Kq \ vc_Vq_outer;
tildeWq = vc_Vq_outer'*vc_Vq_outer;

% if norm(t - tildeWq, 'fro') / norm(t, 'fro') > 1e-6
%   warning('gen_tildeWq: inconsistent t and tildeWq, rel=%.6e', norm(t - tildeWq, 'fro') / norm(t, 'fro'));
% end
% tildeWq = tildeWq + outer_data.tildeVq(:, :, iqibz);



end