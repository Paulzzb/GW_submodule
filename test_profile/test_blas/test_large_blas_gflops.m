% License-Identifier: GPL
%
% Large-matrix BLAS benchmark:
%   1) GEMM: C = A * B, A/B are (n x n)
%   2) TRMM (left): Y = U * X, U upper-triangular (n x n), X (n x rhs)
%   3) TRSM (left): Z = U \ X, U upper-triangular (n x n), X (n x rhs)
%
% Run:
%   cd GW/test_profile/test_blas
%   test_large_blas_gflops

rng(1);

dtype_list = {'double', 'double'};
n_list = [1024, 2048, 4096];
rhs_list = [256, 512, 1024];

% Max timing repeats for stable numbers. Large sizes auto-use fewer repeats.
repeats_small = 3;
repeats_large = 2;

nr = numel(dtype_list) * numel(n_list) * numel(rhs_list);
row = 0;

R = struct( ...
  'dtype', strings(nr, 1), ...
  'n', zeros(nr, 1), ...
  'rhs', zeros(nr, 1), ...
  't_gemm_best_s', zeros(nr, 1), ...
  't_gemm_mean_s', zeros(nr, 1), ...
  'gflops_gemm_best', zeros(nr, 1), ...
  'gflops_gemm_mean', zeros(nr, 1), ...
  't_trmm_best_s', zeros(nr, 1), ...
  't_trmm_mean_s', zeros(nr, 1), ...
  'gflops_trmm_best', zeros(nr, 1), ...
  'gflops_trmm_mean', zeros(nr, 1), ...
  't_trsm_best_s', zeros(nr, 1), ...
  't_trsm_mean_s', zeros(nr, 1), ...
  'gflops_trsm_best', zeros(nr, 1), ...
  'gflops_trsm_mean', zeros(nr, 1));

fprintf('\n=== Large BLAS benchmark (GEMM/TRMM/TRSM) ===\n');
fprintf('dtype = [double, double]\n');
fprintf('n_list = [%s]\n', num2str(n_list));
fprintf('rhs_list = [%s]\n\n', num2str(rhs_list));

for idt = 1:numel(dtype_list)
  dtype = dtype_list{idt};
  fprintf('\n================ dtype = %s ================\n', dtype);

  for in = 1:numel(n_list)
    n = n_list(in);
    rpt = repeats_small;
    if n >= 4096
      rpt = repeats_large;
    end

    fprintf('\n--- n = %d (repeats=%d) ---\n', n, rpt);

    A = rand(n, n, dtype);
    B = rand(n, n, dtype);

    U = triu(rand(n, n, dtype));
    if strcmp(dtype, 'double')
      U = U + eye(n, 'double');
    else
      U = U + eye(n);
    end

    % Warmup GEMM once to reduce first-call effects.
    Cw = A * B; %#ok<NASGU>

    tg = zeros(rpt, 1);
    for ir = 1:rpt
      t0 = tic;
      C = A * B; %#ok<NASGU>
      tg(ir) = toc(t0);
    end
    t_gemm_best = min(tg);
    t_gemm_mean = mean(tg);
    flops_gemm = 2 * n * n * n;
    gflops_gemm_best = flops_gemm / max(t_gemm_best, eps) / 1e9;
    gflops_gemm_mean = flops_gemm / max(t_gemm_mean, eps) / 1e9;

    fprintf('GEMM n=%4d | best %.4fs | mean %.4fs | GFLOPS(best)=%.2f | GFLOPS(mean)=%.2f\n', ...
      n, t_gemm_best, t_gemm_mean, gflops_gemm_best, gflops_gemm_mean);

    for irhs = 1:numel(rhs_list)
      rhs = rhs_list(irhs);
      if rhs > n
        continue;
      end

      row = row + 1;
      X = rand(n, rhs, dtype);

      % Warmup triangular ops.
      Yw = U * X; %#ok<NASGU>
      Zw = U \ X; %#ok<NASGU>

      tt = zeros(rpt, 1);
      for ir = 1:rpt
        t0 = tic;
        Y = U * X; %#ok<NASGU>
        tt(ir) = toc(t0);
      end
      t_trmm_best = min(tt);
      t_trmm_mean = mean(tt);

      ts = zeros(rpt, 1);
      for ir = 1:rpt
        t0 = tic;
        Z = U \ X; %#ok<NASGU>
        ts(ir) = toc(t0);
      end
      t_trsm_best = min(ts);
      t_trsm_mean = mean(ts);

      % Approximate FLOPs for triangular matrix multiply/solve.
      % The exact count differs by implementation details; this is a useful
      % normalization for relative throughput comparison.
      flops_tri = n * n * rhs;
      gflops_trmm_best = flops_tri / max(t_trmm_best, eps) / 1e9;
      gflops_trmm_mean = flops_tri / max(t_trmm_mean, eps) / 1e9;
      gflops_trsm_best = flops_tri / max(t_trsm_best, eps) / 1e9;
      gflops_trsm_mean = flops_tri / max(t_trsm_mean, eps) / 1e9;

      R.dtype(row) = string(dtype);
      R.n(row) = n;
      R.rhs(row) = rhs;
      R.t_gemm_best_s(row) = t_gemm_best;
      R.t_gemm_mean_s(row) = t_gemm_mean;
      R.gflops_gemm_best(row) = gflops_gemm_best;
      R.gflops_gemm_mean(row) = gflops_gemm_mean;
      R.t_trmm_best_s(row) = t_trmm_best;
      R.t_trmm_mean_s(row) = t_trmm_mean;
      R.gflops_trmm_best(row) = gflops_trmm_best;
      R.gflops_trmm_mean(row) = gflops_trmm_mean;
      R.t_trsm_best_s(row) = t_trsm_best;
      R.t_trsm_mean_s(row) = t_trsm_mean;
      R.gflops_trsm_best(row) = gflops_trsm_best;
      R.gflops_trsm_mean(row) = gflops_trsm_mean;

      fprintf('  rhs=%4d | TRMM best %.4fs (%.2f GF/s) | TRSM best %.4fs (%.2f GF/s)\n', ...
        rhs, t_trmm_best, gflops_trmm_best, t_trsm_best, gflops_trsm_best);
    end

    clear A B U Cw;
  end
end

R = trim_results(R, row);
T = table( ...
  R.dtype, R.n, R.rhs, ...
  R.t_gemm_best_s, R.t_gemm_mean_s, R.gflops_gemm_best, R.gflops_gemm_mean, ...
  R.t_trmm_best_s, R.t_trmm_mean_s, R.gflops_trmm_best, R.gflops_trmm_mean, ...
  R.t_trsm_best_s, R.t_trsm_mean_s, R.gflops_trsm_best, R.gflops_trsm_mean, ...
  'VariableNames', { ...
    'dtype', 'n', 'rhs', ...
    't_gemm_best_s', 't_gemm_mean_s', 'gflops_gemm_best', 'gflops_gemm_mean', ...
    't_trmm_best_s', 't_trmm_mean_s', 'gflops_trmm_best', 'gflops_trmm_mean', ...
    't_trsm_best_s', 't_trsm_mean_s', 'gflops_trsm_best', 'gflops_trsm_mean'});

[~, ord_dtype] = ismember(cellstr(T.dtype), dtype_list);
[~, ord] = sortrows([ord_dtype(:), T.n, T.rhs], [1 2 3]);
T = T(ord, :);

fprintf('\n=== Summary ===\n');
disp(T);

function R = trim_results(R, row)
f = fieldnames(R);
for i = 1:numel(f)
  R.(f{i}) = R.(f{i})(1:row, :);
end
end
