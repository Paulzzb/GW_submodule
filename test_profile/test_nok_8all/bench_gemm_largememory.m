function bench_gemm_largememory()
%BENCH_GEMM_LARGEMEMORY Measure double-complex GEMM throughput.
%
% Small dimensions: 128:128:1024
% Larger dimensions: 1x / 4x / 16x of each small size
%
% Run directly in MATLAB:
%   bench_gemm_largememory
%
% Or submit with sbatch:
%   sbatch s_bench_gemm_largememory

  small_sizes = 128:128:1024;
  large_factors = [1, 4, 16];
  n_warmup = 1;
  n_repeat = 2;
  sink = 0.0;

  fprintf('=== GEMM benchmark (double complex) ===\n');
  fprintf('Host: %s\n', char(java.net.InetAddress.getLocalHost.getHostName));
  fprintf('Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
  try
    fprintf('maxNumCompThreads: %d\n', maxNumCompThreads); %#ok<MXMCT>
  catch
    fprintf('maxNumCompThreads: <unavailable>\n');
  end
  fprintf('warmup=%d repeat=%d\n', n_warmup, n_repeat);
  fprintf('FLOP model (complex GEMM): 8*M*N*K\n\n');

  fprintf('%6s %6s %8s %8s %8s %12s %12s %14s %10s\n', ...
    'nSmall', 'xL', 'M', 'K', 'N', 'Best(s)', 'Avg(s)', 'GFLOPS(best)', 'Mem(GB)');
  fprintf('%s\n', repmat('-', 1, 98));

  peak_gflops = 0.0;
  peak_shape = [0, 0, 0];
  peak_desc = '';

  for n_small = small_sizes
    for fac = large_factors
      % Representative shape for "small vs large" dimensions:
      % A: MxK, B: KxN, C: MxN
      m = n_small * fac;
      k = n_small;
      n = n_small * fac;

      % double complex ~= 8 bytes / element.
      mem_gb = 8.0 * (m * k + k * n + m * n) / 1024.0^3;
      A = complex(rand(m, k, 'double'), rand(m, k, 'double'));
      B = complex(rand(k, n, 'double'), rand(k, n, 'double'));

      for iw = 1:n_warmup
        C = A * B; %#ok<NASGU>
      end

      t_all = zeros(n_repeat, 1);
      for ir = 1:n_repeat
        t0 = tic;
        C = A * B;
        t_all(ir) = toc(t0);
        sink = sink + real(C(1, 1));
      end

      t_best = min(t_all);
      t_avg = mean(t_all);
      flops = 8.0 * m * n * k;
      gflops_best = flops / t_best / 1e9;

      if gflops_best > peak_gflops
        peak_gflops = gflops_best;
        peak_shape = [m, k, n];
        peak_desc = sprintf('nSmall=%d, factor=%d', n_small, fac);
      end

      fprintf('%6d %6d %8d %8d %8d %12.6f %12.6f %14.2f %10.2f\n', ...
        n_small, fac, m, k, n, t_best, t_avg, gflops_best, mem_gb);
      clear A B C
    end
  end

  fprintf('\nPeak GFLOPS = %.2f at MxKxN = %dx%dx%d (%s)\n', ...
    peak_gflops, peak_shape(1), peak_shape(2), peak_shape(3), peak_desc);
  fprintf('Checksum (ignore): %.6e\n', sink);
  fprintf('=== end GEMM benchmark ===\n');
end
