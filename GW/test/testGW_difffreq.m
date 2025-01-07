% clear all;
% close all;
dbstop if error
cd ../../../;
addpath(genpath('./'));
randn('state', 0);
rand('state', 0);
KSSOLV_startup;
cd silicon
% si2
si8
cd ../src/GW/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We gives a most simple choice for options_input, check ~/src/ISDF/ISDF.m 
% and ~/src/GW/COmegaCstar.m (or ~/src/GW/README, functions) for details.
options_input.isGW          = true; 
options_input.frequency_dependence      = 2; 
options_input.amin          = 5.0; % Currently, this options_input with unit Ha to consist with the main part of GW method 
options_input.nv            = mol.nel/2;
options_input.nc            = mol.nbnd-options_input.nv-1;
options_input.input         = 'kssolv';
options_input.isbse         = false;
options_input.isISDF        = true;
options_input.iscauchy      = true;
options_input.vcrank_ratio  = 16; 
options_input.vsrank_ratio  = 16; 
options_input.ssrank_ratio  = 16; 
options_input.fileName      = 'GWoutput.mat';
options_input.nv_ener       = mol.nel/2;
options_input.nc_ener       = mol.nel/2;
options_input.nc_ener       = mol.nbnd-options_input.nv-1;
options_input.nv_oper       = mol.nel/2;
options_input.nc_oper       = mol.nel/2;
options_input.nc_oper       = mol.nbnd-options_input.nv-1;

options_input.frequency_dependence_method= 2;
options_input.setting_method = 'kssolv';

options_input.delta_frequency = 2.0;
options_input.frequency_low_cutoff = 40;
options_input.number_imaginary_freqs = 20;   %default number of frequency imaginary part

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following parts are necessary for full frequency approximation
% All are corresponding to inputs in epsilon.inp --> Epsilon/inread.f90
% if options_input.frequency_dependence == 2
% % options_input.delta_frequency_step    = 20.0; % default value of delta frequency, Unit eV.
%   options_input.frequency_low_cutoff = 200; %default value, Unit: eV
%   options_input.broadening = 0.1;
%   if(options_input.frequency_dependence_method == 2)
%     options_input.number_imaginary_freqs = 15;   %default number of frequency imaginary part
%   %  options_input.delta_freq = 20.0;  %default value of delta frequency, Unit: eV
%     options_input.plasma_freq = 2;   %default value of plasma for frequency imaginary part, Unit: Ry
%     options_input.cd_int_method = 0; %For contour deformation calculations, specify the integration method on the imaginary axis (default is 0)
%   elseif(options_input.full_freq_method == 1)
%     error('Error: Not support now!!');
%   elseif(options_input.full_freq_method == 0)
%     options_input.freqevalstep = 0.2; % similar to BGW
%     options_input.freq_cutoff_high = 4.0*options_input.freq_cutoff;
%     options_input.delta_frequency_step = 1.0;
%   else
%     error('Error: Not support now!!');
%   end
  
%   options_input.frequency_dependence_method= 2;
%   options_input.delta_frequency = 20;
% end

% The following parts are necessary for full frequency approximation
% All are corresponding to inputs in epsilon.inp --> Epsilon/inread.f90
% if options_input.frequency_dependence == 2
%   options_input.delta_frequency_eval = 0.2; 
%   options_input.freqevalmin = 0;
%   options_input.number_frequency_eval = 1;
%   options_input.frequency_grid_shift = 2;
%   options_input.max_freq_eval = 2.0;
%   options_input.cd_integral_method = 0;
% end




dbstop if error;
% cd testSi2
cd Si8

% ksinfo.dFreqBrd(1) = 0.0;
options = GWOptions(mol, options_input);
ksinfor = ksinfo(mol, options.Groundstate);
ksinfor_new = setksinfo(ksinfor, options.Groundstate);
% gw_fullfreq_cd(ksinfor, options);
profile on
sintISDF = gw_fullfreq_cd_int_ISDF(ksinfor_new, options)
profile off
saveprof
% sint = gw_fullfreq_cd_int(ksinfor_new, options)
% sres = gw_fullfreq_cd_res(ksinfor_new, options)
% gwCalculation(ksinfo, options;
load GWoutput.mat
% cd ../../../test
disp('testGW done, check in ../src/GW/Si8/GWenergy.mat for GW quasiparticle energies.')
return
