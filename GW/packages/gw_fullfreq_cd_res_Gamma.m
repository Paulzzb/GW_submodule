function sres = gw_fullfreq_cd_res_Gamma(config)
%GW_FULLFREQ_CD_RES_GAMMA  Residual term for single-(k,q) service path.

cleanup = QPlog_push('Fullfreq-CD-Residual-Gamma');

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end

system_data = system.get();
ev = double(system_data.Eo(:, 1, 1)) * ry2ev;
focc = double(system_data.f(:, 1, 1));

nbmin = config.SYSTEM.energy_band_index_min;
nbmax = config.SYSTEM.energy_band_index_max;
bandtocal = nbmin:nbmax;
n_ener = numel(bandtocal);
nb_total = numel(ev);
nsum = min(config.SYSTEM.number_bands_in_summation, nb_total);
nv = find(focc > 1 - TOL_SMALL, 1, 'last');
if isempty(nv)
  error('gw_fullfreq_cd_res_Gamma:occupation', 'Cannot determine nv from occupations.');
end

bandtocal_occ = bandtocal(bandtocal <= nv);
bandtocal_unocc = bandtocal(bandtocal > nv);
nv_ener = numel(bandtocal_occ);

grid_real = config.freqinfo.grid_real;
coeff_real_func = config.freqinfo.coeff_real_func;
nfreq_real = numel(grid_real);

pattern = zeros(n_ener, nsum, nfreq_real);

for ibe = 1:nv_ener
  ibandener = bandtocal_occ(ibe);
  for ibo = 1:nv
    x = (ev(ibandener) - ev(ibo)) + TOL_SMALL;
    if x >= 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      if abs(coeff) > TOL_ZERO
        pattern(ibe, ibo, ifreq) = 1;
      end
    end
  end
end

for ibe = (nv_ener + 1):n_ener
  ibandener = bandtocal_unocc(ibe - nv_ener);
  for ibo = (nv + 1):nsum
    x = (ev(ibandener) - ev(ibo)) + TOL_SMALL;
    if x < 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      if abs(coeff) > TOL_ZERO
        pattern(ibe, ibo, ifreq) = 1;
      end
    end
  end
end

nm_Womega_nm_list = gw_fullfreq_cd_core_Gamma(config, [nbmin, nbmax], [1, nsum], grid_real, pattern);

sres = zeros(n_ener, 1);

occ_sign = -1;
for ibe = 1:nv_ener
  ibandener = bandtocal_occ(ibe);
  for ibo = 1:nv
    x = (ev(ibandener) - ev(ibo)) + TOL_SMALL;
    if x >= 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      if pattern(ibe, ibo, ifreq) > 0
        coeff = coeff_real_func{ifreq}(x);
        sres(ibe) = sres(ibe) - coeff * occ_sign * nm_Womega_nm_list(ibe, ibo, ifreq);
      end
    end
  end
end

occ_sign = 1;
for ibe = (nv_ener + 1):n_ener
  ibandener = bandtocal_unocc(ibe - nv_ener);
  for ibo = (nv + 1):nsum
    x = (ev(ibandener) - ev(ibo)) + TOL_SMALL;
    if x < 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      if pattern(ibe, ibo, ifreq) > 0
        coeff = coeff_real_func{ifreq}(x);
        sres(ibe) = sres(ibe) - coeff * occ_sign * nm_Womega_nm_list(ibe, ibo, ifreq);
      end
    end
  end
end

end
