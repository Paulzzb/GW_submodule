% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function [fftgrid_i, fft_sz, nfft_d, fftgrid_c, Nmu, R_coarse_RLU, R_rot_coarse] = isdf_build_coarse_grid()
% Thin wrapper: use isdf.gen_coeff_coarse('fft_grid') (implementation in +isdf).

  [fftgrid_i, fft_sz, nfft_d, fftgrid_c, Nmu, R_coarse_RLU, R_rot_coarse] = isdf.gen_coeff_coarse('fft_grid');
end
