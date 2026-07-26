% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function print_coarse_grid_report(id)
% Print dense/coarse FFT grid summary (same style as isdftest/gen_indices_coarse_test).
% No-op if pool id is not coarse (interp_scheme ~= "coarse").

  if nargin < 1 || isempty(id)
    id = isdftest.isdftest_current();
  end

  isdf_data = isdftest.get(id);
  if ~strcmp(isdf_data.interp_scheme, "coarse")
    return;
  end

  fft_data = FFT.get();
  fftgrid_i = int32(fft_data.fftgrid(:).');
  fftgrid_c = isdf_data.fftgrid_c;
  Nmu = isdf_data.nisdf;

  N_dense = prod(double(fftgrid_i));
  N_coarse = double(Nmu);

  fprintf('\n=== ISDF coarse-grid pipeline: grid sizes ===\n');
  fprintf('Dense (fine) FFT dims   [%d, %d, %d]\n', fftgrid_i(1), fftgrid_i(2), fftgrid_i(3));
  fprintf('Dense total sites       %d\n', N_dense);
  fprintf('Sparse (coarse) dims    [%d, %d, %d]\n', fftgrid_c(1), fftgrid_c(2), fftgrid_c(3));
  fprintf('Sparse total sites Nmu  %d\n', N_coarse);
  fprintf('Dense / coarse ratio    %.6g\n', N_dense / N_coarse);
  fprintf('Coarse / dense fraction %.6g\n', N_coarse / N_dense);
  fprintf('=============================================\n\n');
end
