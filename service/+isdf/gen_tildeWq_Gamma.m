% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/11 ZZ

function tildeWq = gen_tildeWq_Gamma(id_vc, iqibz, Kq, id_outer, flagherm)
%GEN_TILDEWQ_GAMMA  Build screened ISDF kernel W~ at Gamma.
%
%   tildeWq = isdf.gen_tildeWq_Gamma(id_vc, iqibz, Kq, id_outer)
%   tildeWq = isdf.gen_tildeWq_Gamma(id_vc, iqibz, Kq, id_outer, flagherm)
%
% Forms (this routine computes the K^{-1} contraction; bare V term lives on
% the outer slot from gen_tildeVq when needed by the caller):
%
%   W~(outer,outer) ~ <outer|V_q|vc> K_q^{-1} <vc|V_q|outer>
%
%   id_vc     — ISDF slot used to build K (typically desc 'vc')
%   iqibz     — q IBZ index (Gamma path currently forces iqibz = 1)
%   Kq        — from isdf.gen_Kq_Gamma
%   id_outer  — outer slot: usually 'vn' (SEX); 'nn' for COH / full-freq
%   flagherm  — true: prefer chol(K) (imag-axis / Hermitian K);
%               false: general K\ (real-axis). Default false.

if nargin < 4
  error('gen_tildeWq: Missing input: id_vc, iqibz, Kq, id_outer');
end

if nargin < 5
  flagherm = false;
end

iqibz = 1;

vc_data = isdf.get(id_vc);
outer_data = isdf.get(id_outer);
vc_fac = vc_data.CCHq_trunc_factors{double(iqibz)};
Nkeep_vc = vc_fac.N_keep;
outer_fac = outer_data.CCHq_trunc_factors{double(iqibz)};
Nkeep_outer = outer_fac.N_keep;
coulomb_data = coulomb.get();
vcoul_q = coulomb_data.vcoul(:, iqibz);
if iqibz == 1
  vcoul_q(1) = coulomb_data.vcoul0;
end

tildeWq = zeros(Nkeep_outer, Nkeep_outer);
% <vc|V_q|outer> in G-space (truncated SVD factors)
helperqG_vc = vc_data.helperqG(:, 1:Nkeep_vc, iqibz);
helperqG_outer = outer_data.helperqG(:, 1:Nkeep_outer, iqibz);

vc_Vq_outer = helperqG_vc(:, 1:Nkeep_vc)' * diag(vcoul_q) * helperqG_outer(:, 1:Nkeep_outer);

% Contract with K^{-1}: Hermitian path uses chol when possible
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

end
