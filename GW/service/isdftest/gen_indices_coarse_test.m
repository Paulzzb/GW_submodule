% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function gen_indices_coarse_test()
% GEN_INDICES_COARSE_TEST  ISDF coarse-grid pipeline (driver).
%
% Steps: isdf.coeff.gen_coeff_coarse ('fft_grid', 'wf') -> isdf_build_tildeVq ->
%        isdf_coarse_validate_energies.
%
% Caching of wf_on_coarse / tildeVq is not implemented here (reserved for DB layer).
  fft_data = FFT.get();
  fftgrid_i = int32(fft_data.fftgrid(:).');
  [fftgrid_c, Nmu, R_coarse_RLU, R_rot_coarse] = isdf.coeff.gen_coarse_Rgrid();

  
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

  wf_on_coarse = isdf.coeff.gen_coeff_coarse('wf', fftgrid_c, Nmu, R_coarse_RLU);

  tildeVq = isdf_build_tildeVq(wf_on_coarse, R_coarse_RLU, R_rot_coarse, fftgrid_i, fftgrid_c, fft_sz, Nmu);

  [Esum2, EsumISDF2, DiffEsum2] = isdf_coarse_validate_energies(tildeVq, wf_on_coarse, R_rot_coarse, ...
    R_coarse_RLU, fftgrid_i, fftgrid_c, Nmu);

  fprintf('\n=== ISDF Coarse Validation Report ===\n');
  fprintf('Sum Ex_t^2 (direct):      %.4e\n', Esum2);
  fprintf('Sum Ex_ISDF^2:            %.4e\n', EsumISDF2);
  fprintf('Sum (Ex_t - Ex_ISDF)^2:   %.4e\n', DiffEsum2);
  fprintf('====================================\n\n');
end
