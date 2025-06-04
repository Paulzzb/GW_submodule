function sint = gw_fullfreq_cd_int2(GWinfo, options)
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


startintgral = tic;

% 0. Initialize
% Initialize constant from options.Constant
nameConstants = fieldnames(options.Constant);
for i = 1:numel(nameConstants)
  fieldname = nameConstants{i};
  value = options.Constant.(fieldname);    
  if ~isempty(value)
    strEval = sprintf('%s = %.16f;', fieldname, value);
    eval(strEval);
  end
end
bandtocal = nv-nv_ener+1:nv+nc_ener;
n_ener = length(bandtocal);
ev = GWinfo.ev * ry2ev;

% 1. Generate frequency sequences in imaginary axis.
startFrequencyGen = tic;
[~, ~, grid_imag, coeff_imag_func] = freqgen(GWinfo, options);
nfreq_imag = length(grid_imag);
timeFrequencyGen = toc(startFrequencyGen);
fprintf('Time for Generating Frequencies = %.3s.\n', timeFrequencyGen);




sint = zeros(n_ener, 1);
nm_Womega_nm_list = zeros(n_ener, n_oper, nfreq_imag);

% for ifreq = 1:nfreq_imag
nm_Womega_nm_list(:, :, :) = fourcenterintegral(GWinfo, options, 1, ...
            [nv-nv_ener+1, nv+nc_ener], [nv-nv_oper+1, nv+nc_oper], ...
            grid_imag);
% end




for ibandener = nv-nv_ener+1:nv+nc_ener
  ibe = ibandener-nv+nv_ener;
  for ibandoper = nv-nv_oper+1:nv+nc_oper
    ibo = ibandoper-nv+nv_oper;
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
timeintegral = toc(startintgral);
fprintf("Time to integral on Imaginary axis = %.2d.\n", timeintegral);
end
