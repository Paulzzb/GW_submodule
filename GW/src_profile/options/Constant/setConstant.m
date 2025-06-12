function Constantstruct = setConstant(data, config)

Constantstruct = default_constant();

% Set system based values
nv = data.sys.ne / 2;
nb = length(data.ev);
nc = nb - nv; 
Constantstruct.nv = nv; 
Constantstruct.nc = nc; 
Constantstruct.nv_oper = nv;
Constantstruct.nc_oper = config.SYSTEM.number_bands_in_summation - nv;
Constantstruct.nv_ener = nv - config.SYSTEM.energy_band_index_min + 1;
Constantstruct.nc_ener = config.SYSTEM.energy_band_index_max - nv;

Constantstruct.ne = 'SHIT, ne useful';
Constantstruct.vol = data.sys.vol;
Constantstruct.ng = 'SHIT, nguseful';
Constantstruct.nr = data.sys.n1*data.sys.n2*data.sys.n3;

end % EOF