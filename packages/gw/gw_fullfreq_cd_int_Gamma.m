% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/13 ZZ

function sint = gw_fullfreq_cd_int_Gamma(config)
%GW_FULLFREQ_CD_INT_GAMMA  Imaginary-axis integral term for double-(k,q) service path.

cleanup = output.push('Fullfreq-CD-Integral-Gamma');

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end

system_data = system.get();
ev = double(system_data.Eo(:, 1, 1)) * ry2ev;

nbmin = config.SYSTEM.energy_band_index_min;
nbmax = config.SYSTEM.energy_band_index_max;
bandtocal = nbmin:nbmax;
n_ener = numel(bandtocal);
nb_total = numel(ev);
nsum = min(config.SYSTEM.number_bands_in_summation, nb_total);

nfreq_imag = config.FREQUENCY.number_imaginary_freqs;
grid_imag = config.freqinfo.grid_imag;
coeff_imag_func = config.freqinfo.coeff_imag_func;

pattern = ones(n_ener, nsum, nfreq_imag);
nm_Womega_nm_list = gw_fullfreq_cd_core_Gamma(config, [nbmin, nbmax], [1, nsum], grid_imag, pattern);

sint = zeros(n_ener, 1);
for ibe = 1:n_ener
  ibandener = bandtocal(ibe);
  for ibo = 1:nsum
    x = (ev(ibandener) - ev(ibo)) + TOL_SMALL;
    for ifreq = 1:nfreq_imag
      coeff = coeff_imag_func{ifreq}(x);
      sint(ibe) = sint(ibe) + coeff * nm_Womega_nm_list(ibe, ibo, ifreq);
    end
  end
end

sint = sint / pi;

end
