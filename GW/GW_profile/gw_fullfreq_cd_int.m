function sint = gw_fullfreq_cd_int(GWinfo, config)
% This file will be final version for the full-frequency GW calculation.
% Generally speaking, user should not directly call this function, instead
% call this indirectly by calling through gwCalculation->gw_fullfreq_cd.

% Description
%     Calculate the self-energies using contour integral idea within the GW method.
%     This code only computes the integral on the imaginary axis,  

% Parameters
%   Input:
%     GWinfo: class @GWinfo, contains ground-state information.
%     options: class @GWOptions, contains necessary parameters for the calculation.
%   Output:
%     sint: the integral on imaginary axis, size nbandener*1.

% Structure 
%   The code is organized as followed:
%     0. Initialize inputs.
%     1. Generate the qudrature weights and the quadrature node on the nary axis.
%        For each frequencies omega
%     2    Calculate W(omega)
%          For each band to calculate n
%     3.     Calculate <nn'|W(omega)|nn'> for all n'. 
%     4.     Multiply the coefficient and add the result to sint  
%          end 
%        end 
cleanup = QPlog_push('Fullfreq-CD-Integral');


startintgral = tic;

% 0. Initialize
default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end

nbmin = config.SYSTEM.energy_band_index_min;
nbmax = config.SYSTEM.energy_band_index_max;
bandtocal = nbmin:nbmax;
n_ener = length(bandtocal);
n_oper = config.SYSTEM.number_bands_in_summation;
nsum = config.SYSTEM.number_bands_in_summation;

ev = GWinfo.ev * ry2ev;

nfreq_imag = config.FREQUENCY.number_imaginary_freqs;
grid_imag = config.freqinfo.grid_imag;
coeff_imag_func = config.freqinfo.coeff_imag_func;





sint = zeros(n_ener, 1);
nm_Womega_nm_list = zeros(n_ener, n_oper, nfreq_imag);

nm_Womega_nm_list(:, :, :) = fourcenterintegral(GWinfo, config, 1, ...
            [nbmin, nbmax], [1, nsum], grid_imag);

for ibandener = nbmin:nbmax
  ibe = ibandener-nbmin+1;
  for ibandoper = 1:nsum 
    ibo = ibandoper;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    for ifreq = 1:nfreq_imag
      coeff_func = coeff_imag_func{ifreq};
      coeff = coeff_func(x);
      nm_W_nm = nm_Womega_nm_list(ibe, ibo, ifreq);
      sint(ibe) = sint(ibe) + coeff*nm_W_nm;
    end
  end
end

% sint = sint / pi * ry2ev;
sint = sint / pi;
% timeintegral = toc(startintgral);
% sprintf("Time to integral on Imaginary axis = %.2d.\n", timeintegral);
end
