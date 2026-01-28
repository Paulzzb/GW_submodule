% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/28 ZZ
%
function out = isdf_desc(GWinfo, config)
% Create describer for current input
  out = struct();
  out.desc_sys    = desc();
  out.desc_param  = desc();
  out.desc_type1  = desc();
  out.desc_type2  = desc();
  out.desc_type3  = desc();
  %
  out.schema = struct();
  out.schema.key_fields = struct( ...
    'desc_sys',   true, ...
    'desc_param', true, ...
    'desc_type1', true, ...
    'desc_type2', true, ...
    'desc_type3', true ...
  ); 
  % -----------------------------------------------------------------
  desc_sys.add("prefix", config.CONTROL.prefix);
  desc_sys.add("coulomb_truncation_method", config.CUTOFFS.coulomb_truncation_method);
  desc_sys.add("coulomb_truncation_parameter", config.CUTOFFS.coulomb_truncation_parameter);
  desc_sys.add("coulomb_cutoff", config.CUTOFFS.coulomb_cutoff);
  bmax1 = config.SYSTEM.number_bands_max;
  bmin1 = config.SYSTEM.number_bands_min;
  bmax2 = config.SYSTEM.energy_band_index_max;
  bmin = config.SYSTEM.energy_band_index_min;
  bmax = max(bmax1, bmax2);
  bmin = min(bmin1, bmin);
  desc_sys.add("bands_max", bmax);
  desc_sys.add("bands_min", bmin);
  desc_sys.add("nr", config.ISDF.sys.nr);
  desc_sys.add("ng", config.ISDF.sys.ng);
  desc_sys.add("vol", config.ISDF.sys.vol);
  desc_sys.add("supercell", config.ISDF.sys.supercell);
  % In order to avoid that, another groundstate calculation is performed, but
  % forget to remove the old ISDF database, we need some comparison of GS.
  % Of course I cannot save all groundstate info for compare, so I save energy
  % structure only.
  desc_sys.add("groundstate_energy", GWinfo.ev);
  %
  out.desc_sys = desc_sys;
  % -----------------------------------------------------------------
  desc_param.add("exxmethod", config.ISDF.exxmethod);
  desc_param.add("seed", config.ISDF.seed);
  desc_param.add("init", config.ISDF.init);
  desc_param.add("weight", config.ISDF.weight);
  out.desc_param = desc_param;
  % -----------------------------------------------------------------
  descA = desc();
  descA.add('isdf_ratio_type1', config.ISDF.isdf_ratio_type1);
  out.descA = descA;
  % -----------------------------------------------------------------
  descB = desc();
  descB.add('isdf_ratio_type2', config.ISDF.isdf_ratio_type2);
  out.descB = descB;
  % -----------------------------------------------------------------
  descC = desc();
  descC.add('isdf_ratio_type3', config.ISDF.isdf_ratio_type3);
  out.descC = descC;
end
