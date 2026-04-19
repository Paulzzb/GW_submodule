% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function FFT_G_table()
  %
  FFT.getnGo()
  %
  r_lat_data = lattice.manager('r_lat', 'get');
  Ggrid = r_lat_data.Ggrid_RLU;
  %
  FFT_data = FFT.get();
  fftgrid = FFT_data.fftgrid;
  n1 = fftgrid(1);
  n2 = fftgrid(2); 
  n3 = fftgrid(3);
  %
  for iGo = 1:FFT_data.nGo
    G_Go = Ggrid - Ggrid(iGo, :);
    ind = mod(G_Go(:, 1), n1) + mod(G_Go(:, 2), n2) * n1 ...
          + mod(G_Go(:, 3), n3) * n1 * n2 + 1;
    tmp = int32(ind);
    FFT_data.G_table(:, iGo) = tmp;
  end

  FFT.save2mod(FFT_data);
end
