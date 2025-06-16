function sres = gw_fullfreq_cd_res(GWinfo, config)

% This file will be final version for the full-frequency GW calculation.
% Generally speaking, user should not directly call this function, instead
% call this indirectly by calling through gwCalculation->gw_fullfreq_cd.

% Description
%     Calculate the self-energies using contour integral idea within the GW method.
%     This code only computes the summation of the residual.

% Parameters
%   Input:
%     GWinfo: class @GWinfo, contains ground-state information.
%     options: class @GWOptions, contains necessary parameters for the calculation.
%   Output:
%     sres: the summation of the residual, size nbandener*1.

% Structure 
%   The code is organized as followed:
%     0. Initialize inputs.
%     1. Generate the qudrature weights and the quadrature node on the real axis.
%        For each frequencies omega
%     2. Calculate pattern, i.e., for which set of (n, n') we need to 
%        calculate <nn'|W(omega)|nn'>.
%     3    Calculate W(omega)
%     4.   Calculate <nn'|W(omega)|nn'> for all n and n' both occupied
%          and both unoccupied, then save it. 
%          For each n and n' both occupied and unoccupied
%     5.     Check if it is residual, add to sres.
%          end 
%        end 

% 0. Initialize
% Initialize constant from options.Constant
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
nv = find(GWinfo.occupation > 1 - TOL_SMALL, 1, 'last');
if (nv >= nsum)
  msg = sprintf('Number of valence bands = %d >= Number of bands in summation = %d', nv, nsum);
  GWerror(msg)
end


bandtocal_occ = bandtocal(find(bandtocal <= nv));
bandtocal_unocc = bandtocal(find(bandtocal > nv));
nv_ener = length(bandtocal_occ);
nc_ener = length(bandtocal_unocc);

% Initialize other values
% We use unit as ev, while usual input is Ry.
ev = GWinfo.ev * ry2ev;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid_real = config.freqinfo.grid_real;
coeff_real_func = config.freqinfo.coeff_real_func;
nfreq_real = length(grid_real);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nm_Womega_nm_list = zeros(n_ener, n_oper, nfreq_real);
nm_Womega_nm_pattern = zeros(n_ener, n_oper, nfreq_real);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate pattern first, deciding which part to calculate
% After that, nm_Woemga_nm_pattern contains 1 or 0.
% We only need values for place 1.
for ibe = 1:nv_ener 
  ibandener = bandtocal_occ(ibe);
  for ibo = 1:nv
    ibandoper = ibo;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    if x >= 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      if abs(coeff) > TOL_ZERO
        nm_Womega_nm_pattern(ibe, ibo, ifreq) = 1;
      end
    end
  end
end

% Then both unoccupied states
for ibe = nv_ener+1:n_ener
  ibandener = bandtocal_unocc(ibe-nv_ener);
  for ibo = nv+1:nsum
    ibandoper = ibo;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    if x < 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      if abs(coeff) > TOL_ZERO
        nm_Womega_nm_pattern(ibe, ibo, ifreq) = 1;
      end
    end 
  end
end

% Pattern calculated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate four-center integral
nm_Womega_nm_list(:, :, :) = fourcenterintegral(GWinfo, config, 1, ...
            [nbmin, nbmax], [1, nsum], grid_real, ...
            'pattern', nm_Womega_nm_pattern);

% Four-center integral calculated 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate residual energies 
sres = zeros(n_ener, 1);
% Use the elements in MWM to calculate sres
% First, both occupied states.

occ_sign = -1;
for ibe = 1:nv_ener 
  ibandener = bandtocal_occ(ibe);
  for ibo = 1:nv
    ibandoper = ibo;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    if x >= 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      if (nm_Womega_nm_pattern(ibe, ibo, ifreq) > 0)
        coeff = coeff_real_func{ifreq}(x);
        sres(ibe) = sres(ibe) ...
        - coeff * occ_sign * nm_Womega_nm_list(ibe, ibo, ifreq); 
      end
    end
  end
end

% Then both unoccupied states
occ_sign = 1;
for ibe = nv_ener+1:n_ener 
  ibandener = bandtocal_unocc(ibe-nv_ener);
  for ibo = nv+1:nsum
    ibandoper = ibo;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    if x < 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      if (nm_Womega_nm_pattern(ibe, ibo, ifreq) > 0)
        coeff = coeff_real_func{ifreq}(x);
        sres(ibe) = sres(ibe) ...
        - coeff * occ_sign * nm_Womega_nm_list(ibe, ibo, ifreq); 
      end
    end
  end
end

end % function
