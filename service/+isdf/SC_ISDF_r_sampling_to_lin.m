% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/11

function lin = SC_ISDF_r_sampling_to_lin(isdf_data)
%SC_ISDF_R_SAMPLING_TO_LIN  Map R_sampling_RLU rows to fine FFT linear indices.

  fft_data = FFT.get();
  Nmu = double(isdf_data.nisdf);
  g = double(fft_data.fftgrid(:)).';
  Rgrid = double(fft_data.Rgrid_RLU);
  lin = zeros(Nmu, 1, 'int32');
  Rs = double(isdf_data.R_sampling_RLU(1:Nmu, :));
  for i = 1:Nmu
    v = mod(round(Rs(i, :)), g);
    [is_hit, k] = ismember(v, Rgrid, 'rows');
    if ~is_hit
      dd = zeros(size(Rgrid, 1), 3);
      for ddim = 1:3
        t = abs(Rgrid(:, ddim) - v(ddim));
        dd(:, ddim) = min(t, min(abs(t - g(ddim)), abs(t + g(ddim))));
      end
      [~, k] = min(sum(dd, 2));
    end
    lin(i) = int32(k);
  end
end
