% License-Identifier: GPL
%
% Simple benchmark for isdf.prod-like kernel:
%   out = (conj(Psi) * conj(psi)') .* (Phi * phi')
%
% Compare batch sizes used in adaptive_weight:
%   batch = 16, 64, 128, 256
%
% Run:
%   cd GW/test_profile/test_blas
%   test_isdf_prod_batch_gflops

rng(1);

% "Nold" is the second matrix row count in adaptive_weight (Psi_on_grid(1:Nold,:)).
nold_list = [1e4, 5e4, 1e5, 2e5];
batch_list = [16, 64, 128, 256];

% Typical narrow widths (similar to band subspace sizes in isdf workflows).
k1 = 64;
k2 = 64;

% Keep memory moderate for large Nold.
dtype = 'single';  % switch to 'double' if needed
repeats = 3;

nr = numel(nold_list) * numel(batch_list);
row = 0;

results = struct( ...
  'Nold', zeros(nr, 1), ...
  'batch', zeros(nr, 1), ...
  'k1', zeros(nr, 1), ...
  'k2', zeros(nr, 1), ...
  't_best_s', zeros(nr, 1), ...
  't_mean_s', zeros(nr, 1), ...
  'gflops_best', zeros(nr, 1), ...
  'gflops_mean', zeros(nr, 1));

fprintf('\n=== isdf.prod-like kernel benchmark ===\n');
fprintf('dtype=%s, k1=%d, k2=%d, repeats=%d\n', dtype, k1, k2, repeats);
fprintf('kernel: (conj(Psi)*conj(psi)'') .* (Phi*phi'')\n\n');

for in = 1:numel(nold_list)
  Nold = nold_list(in);
  fprintf('--- Nold = %d ---\n', Nold);
  for ib = 1:numel(batch_list)
    bsz = batch_list(ib);
    row = row + 1;

    Psi = rand(bsz, k1, dtype);
    psi = rand(Nold, k1, dtype);
    Phi = rand(bsz, k2, dtype);
    phi = rand(Nold, k2, dtype);

    % warmup
    out_w = (conj(Psi) * conj(psi')) .* (Phi * phi'); %#ok<NASGU>

    tt = zeros(repeats, 1);
    for ir = 1:repeats
      t0 = tic;
      out = (conj(Psi) * conj(psi')) .* (Phi * phi'); %#ok<NASGU>
      tt(ir) = toc(t0);
    end

    t_best = min(tt);
    t_mean = mean(tt);

    % FLOPs estimate:
    % GEMM1: (b x k1) * (k1 x Nold) -> ~2*b*k1*Nold
    % GEMM2: (b x k2) * (k2 x Nold) -> ~2*b*k2*Nold
    % Hadamard: b*Nold
    flops = 2 * bsz * k1 * Nold + 2 * bsz * k2 * Nold + bsz * Nold;

    gflops_best = flops / max(t_best, eps) / 1e9;
    gflops_mean = flops / max(t_mean, eps) / 1e9;

    results.Nold(row) = Nold;
    results.batch(row) = bsz;
    results.k1(row) = k1;
    results.k2(row) = k2;
    results.t_best_s(row) = t_best;
    results.t_mean_s(row) = t_mean;
    results.gflops_best(row) = gflops_best;
    results.gflops_mean(row) = gflops_mean;

    fprintf('batch=%3d | best %.4fs | mean %.4fs | GFLOPS(best)=%.2f | GFLOPS(mean)=%.2f\n', ...
      bsz, t_best, t_mean, gflops_best, gflops_mean);
  end
  fprintf('\n');
end

T = table(results.Nold, results.batch, results.k1, results.k2, ...
  results.t_best_s, results.t_mean_s, results.gflops_best, results.gflops_mean, ...
  'VariableNames', {'Nold', 'batch', 'k1', 'k2', ...
                    't_best_s', 't_mean_s', 'gflops_best', 'gflops_mean'});

[~, ord] = sortrows([T.Nold, T.batch], [1, 2]);
T = T(ord, :);

fprintf('=== Summary ===\n');
disp(T);

