% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/11 ZZ

function report = verify_ellipk_replacement()
%VERIFY_ELLIPK_REPLACEMENT  Compare Driscoll ellipkkp vs MATLAB ellipke.
%
%   report = isdf.Cauchy.verify_ellipk_replacement()
%
% COmegaCstar uses:
%   L = -log(k)/pi;  [K,~] = ellipkkp(L)
% which matches MATLAB built-in K = ellipke(k^2) for K (not Kp).

  ratios = [1.01, 1.1, 1.5, 2, 5, 10, 50, 100, 1e3, 1e4, 1e6, 1e8];
  k_from_ratio = (sqrt(ratios) - 1) ./ (sqrt(ratios) + 1);
  k_grid = unique([ ...
    1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, ...
    0.8, 0.9, 0.95, 0.99, 0.999, 0.9999, 1 - 1e-8, ...
    k_from_ratio], 'stable');

  n = numel(k_grid);
  rel = nan(n, 1);
  rows = cell(n, 1);

  for i = 1:n
    k = k_grid(i);
    if ~(k > 0 && k < 1)
      continue
    end
    L = -log(k) / pi;
    [K_kkp, ~] = isdf.Cauchy.ellipkkp(L);
    [K_ke, ~] = ellipke(k^2);
    rel(i) = abs(K_kkp - K_ke) / max(abs(K_ke), eps);
    rows{i} = sprintf('k=%.6g  L=%.6g  |dK/K|=%.3e', k, L, rel(i));
  end

  mask = isfinite(rel);
  report = struct();
  report.k = k_grid(:);
  report.rel_K = rel;
  report.max_rel_K = max(rel(mask));
  report.ok = report.max_rel_K < 1e-12;
  report.lines = rows(mask);

  fprintf('\n=== verify_ellipk_replacement (ellipkkp vs ellipke) ===\n');
  fprintf('max |dK/K|: %.3e  (ok<%.0e: %d)\n\n', report.max_rel_K, 1e-12, report.ok);
  for i = 1:numel(report.lines)
    fprintf('%s\n', report.lines{i});
  end
  fprintf('\n');
end
