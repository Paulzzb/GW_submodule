function [chi0] = get_chi0(gme_q)
% 初始化累加器
chi0_sum = 0;

% 将所有非空的 gme_q{i} 组合成一个矩阵
X = [gme_q{~cellfun(@isempty, gme_q)}];
X = X(:, any(X ~= 0));  % 去除全零列

% 一次性计算 chi0_sum
chi0_sum = conj(X) * X.';

% 返回结果（带负号）
chi0 = -chi0_sum;
fprintf('Final chi0 calculation completed.\n');
end