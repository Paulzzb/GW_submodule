% 这个脚本用于对于新设计算法的测试工作
%
% 算法说明（TeX）: demo_isdf_algorithm.tex
%
% 实现位置:
%   +isdf   — 现有算法 (helperqR = MCHq / CCHq,  应用: c' * tildeVq * c)
%   +isdftest — SVD 测试算法 (helperqR = MCHq * C^{-1/2};
%               get_rho_xalpha 内做 c_t = C^{-1/2} c，再 c_t' * tildeVq * c_t)
%
% 典型调用流程 (+isdftest):
%   id = isdftest.isdftest_add('test');
%   isdftest.set_nrange(id, config.SYSTEM);
%   isdftest.coeff.gen_coeff(id);
%   isdftest.gen_tildeVq(id);              % 可选第二参 s_cut
%   c_t = isdftest.get_rho_xalpha(id, param);   % 已含 C^{-1/2} 变换
%   E_new = real(c_t' * isdftest.get(id).tildeVq(:,:,iqibz) * c_t);
%   E_old = real(c_rho' * isdf.get(id_vc).tildeVq(:,:,iqibz) * c_rho);  % c_rho 来自 isdf.get_rho_xalpha
%
% 验证目标:
%   1. pseudoinverse 截断阈值 s 对 helper / tildeVq 的影响 (调节 gen_tildeVq 的 s_cut)
%   2. 新 SVD 算法与旧 C^{-1} 路径在能量收缩上是否一致
