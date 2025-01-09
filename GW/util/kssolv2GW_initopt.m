function [options_in, sys] = kssolv2GW_initopt(options_in, mol)

% Initialize options_in with KSSOLV inputs
% Input:
%   mol: KSSOLV class @Moleculer, contains molecule information
%   options_in: Parameters setting for GW calculation given by the users.
% Output:
%   options_in: Parameters setting for GW calculation.
%               If required parameters are not specified by the user,
%               default values will be used.
%   sys: system information, contains
%        ng, nr, ne, vol, xyzlist, n1, n2, n3, supercell.
  
  F = KSFFT(mol);
  [sys.ng, sys.nr] = size(F);
  clear F
  sys.vol = mol.vol;
  sys.xyzlist = mol.xyzlist;
  sys.n1 = mol.n1;
  sys.n2 = mol.n2;
  sys.n3 = mol.n3;
  sys.ne = mol.nel;
  sys.supercell = mol.supercell;

  if ~isfield(options_in, 'isGW')
    options_in.isGW = true;
  end
  if ~isfield(options_in, 'isBSE')
    options_in.isBSE = false;
  end
  if ~isfield(options_in, 'isISDF')
    options_in.isISDF = false;
  end
  if ~isfield(options_in, 'isCauchy')
    if (options_in.isISDF == true)
      options_in.isCauchy = true;
    else
      options_in.isCauchy = false;
    end
  end
  
  
  if ~isfield(options_in, 'iskpoints')
    options_in.iskpoints = false;
  end
  if ~isfield(options_in, 'Spin_dependence')
    options_in.Spin_dependence = 0;
  end
  if ~isfield(options_in, 'frequency_dependence')
    options_in.frequency_dependence = 0;
  end

  if ~isfield(options_in, 'nv')
    if options_in.Spin_dependence == 0
      options_in.nv = mol.nel / 2;
    else
      options_in.nv = 0;
      error('Only support spin degeneration calculation now!');
    end
  end
  if ~isfield(options_in, 'amin')
    options_in.amin = 5; 
  end
  if ~isfield(options_in, 'coulomb_truncation')
    options_in.coulomb_truncation = 2;
  end
  if ~isfield(options_in, 'nc')
    options_in.nc = options_in.nv;
  end
  if ~isfield(options_in, 'input')
    options_in.input = 'kssolv';
  end
  if ~isfield(options_in, 'setting_method')
    options_in.setting_method = options_in.input;
  end
  if ~isfield(options_in, 'dir')
    options_in.dir = [];
  end
  if ~isfield(options_in, 'nv_ener')
    options_in.nv_ener = options_in.nv;
  end
  if ~isfield(options_in, 'nc_ener')
    options_in.nc_ener = options_in.nc;
  end
  if ~isfield(options_in, 'nv_oper')
    options_in.nv_oper = options_in.nv;
  end
  if ~isfield(options_in, 'nc_oper')
    options_in.nc_oper = options_in.nc;
  end
  options_in.n_ener = options_in.nv_ener + options_in.nc_ener;
  options_in.n_oper = options_in.nv_oper + options_in.nc_oper;

  % nv_ener and nv_oper should be less than nv,
  % nc_ener and nc_oper should be less than nc.
  if options_in.nv_ener > options_in.nv
    error("nv_ener should be less than nv");
  end
  if options_in.nc_ener > options_in.nc
    error("nc_ener should be less than nc");
  end
  if options_in.nv_oper > options_in.nv
    error("nv_oper should be less than nv");
  end
  if options_in.nc_oper > options_in.nc
    error("nc_oper should be less than nc");
  end
  
  % GPP initializing
  % if (options_in.frequency_dependence == 1)
% 
  % end
   
  
  % full-frequency initializing
  if (options_in.frequency_dependence == 2)
    if ~isfield(options_in, 'delta_frequency_step')
      options_in.delta_frequency_step = 1; % Need change before release.
    end
    if ~isfield(options_in, 'frequency_low_cutoff')
      options_in.frequency_low_cutoff = 200;
    end
    if ~isfield(options_in, 'imag_freq')
       options_in.imag_freq = 2.0;
    end
    if ~isfield(options_in, 'number_imaginary_freqs')
      options_in.number_imaginary_freqs = 15;
    end
    if ~isfield(options_in, 'broadening')
      options_in.broadening = 0.2;
    end
    if ~isfield(options_in, 'delta_frequency')
      options_in.delta_frequency = 1;
    end
  end




  % if ~isfield(options_in. 'nv_oper')
    % nv_oper = nv;
  % else
    % nv_oper = options.nv_oper;
  % end
  % if ~isfield(options, 'nc_oper')
    % nc_oper = nc;
  % else
    % nc_oper = options.nc_oper;
  % end
 

end
