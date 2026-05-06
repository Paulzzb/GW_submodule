% License-Identifier: GPL
%
% Narrow-BLAS micro benchmark for current server.
% Tests:
%   1) GEMM (A' * B), where A/B are (N x w)
%   2) Triangular multiply (U * B), U upper-triangular (w x w), B (w x N)
%   3) Triangular solve left (U \ B), B (w x N)
%   4) Triangular solve right (B / U), B (N x w)
%
% "Narrow" widths: w = [8, 16, 32, 64]
% Large dimension: N = [1e4, 1e5, 5e5, 1e6]
% dtypes:
%   - double (real)
%   - single (real)
%   - cdouble (complex double)
%   - csingle (complex single)
%
% Run:
%   cd GW/test_profile/test_blas
%   test_narrow_blas

rng(1);

width_list = [8, 16, 32, 64];
n_list = [1e4, 1e5, 5e5, 1e6];
dtype_list = {'double', 'single', 'cdouble', 'csingle'};

nr = numel(width_list) * numel(n_list) * numel(dtype_list);
row = 0;
results = struct( ...
  'dtype', strings(nr, 1), ...
  'w', zeros(nr, 1), ...
  'N', zeros(nr, 1), ...
  't_gemm', zeros(nr, 1), ...
  't_trmm', zeros(nr, 1), ...
  't_trsm_l', zeros(nr, 1), ...
  't_trsm_r', zeros(nr, 1), ...
  'gflops_gemm', zeros(nr, 1), ...
  'gflops_trmm', zeros(nr, 1), ...
  'gflops_trsm_l', zeros(nr, 1), ...
  'gflops_trsm_r', zeros(nr, 1));

fprintf('\n=== Narrow BLAS benchmark ===\n');
fprintf('dtype_list = [double, single, cdouble, csingle]\n');
fprintf('width_list = [%s]\n', num2str(width_list));
fprintf('N_list = [%s]\n\n', num2str(n_list));

for in = 1:numel(n_list)
  N = n_list(in);
  fprintf('\n================ N = %d ================\n', N);
  for idt = 1:numel(dtype_list)
    dtype = dtype_list{idt};
    fprintf('\n--- dtype = %s ---\n', dtype);
    for iw = 1:numel(width_list)
      w = width_list(iw);
      row = row + 1;

      U = make_rand(w, w, dtype);
      U = triu(U);
      if dtype(1) == 'c'
        U = U + eye(w) + 1i * eye(w);
      elseif strcmp(dtype, 'single')
        U = U + eye(w, 'single');
      else
        U = U + eye(w);
      end

      A = make_rand(N, w, dtype);
      B = make_rand(N, w, dtype);
      BL = make_rand(w, N, dtype);
      BR = make_rand(N, w, dtype);

      % warmup
      Cw = A' * B; %#ok<NASGU>
      Yw = U * BL; %#ok<NASGU>
      Xlw = U \ BL; %#ok<NASGU>
      Xrw = BR / U; %#ok<NASGU>

      t0 = tic;
      C = A' * B; %#ok<NASGU>
      t_gemm = toc(t0);

      t0 = tic;
      Y = U * BL; %#ok<NASGU>
      t_trmm = toc(t0);

      t0 = tic;
      Xl = U \ BL; %#ok<NASGU>
      t_trsm_l = toc(t0);

      t0 = tic;
      Xr = BR / U; %#ok<NASGU>
      t_trsm_r = toc(t0);

      flops_gemm = 2 * N * w * w;
      flops_tri = N * w * w;

      results.dtype(row) = string(dtype);
      results.w(row) = w;
      results.N(row) = N;
      results.t_gemm(row) = t_gemm;
      results.t_trmm(row) = t_trmm;
      results.t_trsm_l(row) = t_trsm_l;
      results.t_trsm_r(row) = t_trsm_r;
      results.gflops_gemm(row) = flops_gemm / max(t_gemm, eps) / 1e9;
      results.gflops_trmm(row) = flops_tri / max(t_trmm, eps) / 1e9;
      results.gflops_trsm_l(row) = flops_tri / max(t_trsm_l, eps) / 1e9;
      results.gflops_trsm_r(row) = flops_tri / max(t_trsm_r, eps) / 1e9;

      fprintf('w=%3d | GEMM %.4fs | TRMM %.4fs | TRSM(L) %.4fs | TRSM(R) %.4fs\n', ...
        w, t_gemm, t_trmm, t_trsm_l, t_trsm_r);

      clear A B BL BR C Y Xl Xr Cw Yw Xlw Xrw;
    end
  end
end

T = table( ...
  results.dtype, results.w, results.N, ...
  results.t_gemm, results.gflops_gemm, ...
  results.t_trmm, results.gflops_trmm, ...
  results.t_trsm_l, results.gflops_trsm_l, ...
  results.t_trsm_r, results.gflops_trsm_r, ...
  'VariableNames', { ...
  'dtype', 'w', 'N', ...
  't_gemm_s', 'gflops_gemm', ...
  't_trmm_s', 'gflops_trmm', ...
  't_trsm_l_s', 'gflops_trsm_l', ...
  't_trsm_r_s', 'gflops_trsm_r'});

fprintf('\n=== Summary table ===\n');
% Keep output grouped by N -> dtype -> w
[~, ord_dtype] = ismember(cellstr(T.dtype), dtype_list);
[~, ord] = sortrows([T.N, ord_dtype(:), T.w], [1 2 3]);
T = T(ord, :);
disp(T);

function X = make_rand(m, n, dtype)
switch dtype
  case 'double'
    X = rand(m, n);
  case 'single'
    X = rand(m, n, 'single');
  case 'cdouble'
    X = rand(m, n) + 1i * rand(m, n);
  case 'csingle'
    X = complex(rand(m, n, 'single'), rand(m, n, 'single'));
  otherwise
    error('Unknown dtype: %s', dtype);
end
end

