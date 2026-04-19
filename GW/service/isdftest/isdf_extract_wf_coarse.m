% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function wf_on_coarse = isdf_extract_wf_coarse(fftgrid_i, fft_sz, fftgrid_c, Nmu, R_coarse_RLU)
% Thin wrapper: isdf.gen_coeff_coarse('wf', ...).

  wf_on_coarse = isdf.gen_coeff_coarse('wf', fftgrid_i, fft_sz, fftgrid_c, Nmu, R_coarse_RLU);
end
