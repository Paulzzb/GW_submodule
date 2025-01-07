function sres = gw_fullfreq_cd_res(GWinfo, options)

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
%     2    Calculate W(omega)
%     3.   Calculate <nn'|W(omega)|nn'> for all n and n' both occupied
%          and both unoccupied, then save it. 
%          For each n and n' both occupied and unoccupied
%     4.     Check if it is residual, add to sres.
%          end 
%        end 


testflag1 = true;
testflag1 = false;

startintgral = tic;

flagry2ev = true;
% flagry2ev = false;

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
mi = sqrt(-1);
bandtocal = nv-nv_ener+1:nv+nc_ener;
bandtocal_occ = find(bandtocal <= nv);
bandtocal_unocc = find(bandtocal > nv);
n_ener = length(bandtocal);
nv_ener = length(bandtocal_occ);
nc_ener = length(bandtocal_unocc);

% Initialize other values
if flagry2ev
  GWinfo.ev = GWinfo.ev * ry2ev;
  GWinfo.Vxc = GWinfo.Vxc * ry2ev;
  options.GWCal.dFreqBrd = options.GWCal.dFreqBrd * ry2ev;
  options.GWCal.dFreqGrid = options.GWCal.dFreqGrid * ry2ev;
end

GWinfo.Z     = GWinfo.Z * sqrt(vol);
Z     = GWinfo.Z;
ev    = GWinfo.ev;
Vxc   = GWinfo.Vxc;
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
aqsFlag = ~isempty(GWinfo.aqs);
Dcoul(1, 1) = GWinfo.coulG0;
if flagry2ev
  Dcoul = Dcoul * ry2ev;
end


% 1. Generate frequency sequences in real axis and imaginary axis.
% method_real = 
% method_imag = 
startFrequencyGen = tic;
[grid_real, coeff_real_func, grid_imag, coeff_imag_func] = freqgen(GWinfo, options);
nfreq_real = length(grid_real);
nfreq_imag = length(grid_imag);
% DOUBT!!!
% eta = options.GWCal.dBrdning / ry2ev;
eta = 0.0;
timeFrequencyGen = toc(startFrequencyGen);
fprintf('Time for Generating Frequencies = %.3s.\n', timeFrequencyGen);

% Energies use unit 'ev' in this code.
% Since both KSSOLV_dft and frequency generating part use Ry as unit
% Change unit first.
% ev = ev * ry2ev;
% grid_real = grid_real * ry2ev;
% grid_imag = grid_imag * ry2ev;
% Dcoul = Dcoul * ry2ev;


sres = zeros(n_ener, 1);
if testflag1
  load I_eps_array.mat
  load coeff_real.mat
  load sres_elements.mat  
end

% Place to save the <nn'|W(omega)|nn'>, MWM_occ for two occupied,
% and MWM_unocc for two unoccupied.
MWM_occ = zeros(nv_ener, nv_oper, nfreq_real);
MWM_unocc = zeros(nc_ener, nc_oper, nfreq_real);

for ifreq = 1:nfreq_real
  omega = grid_real(ifreq);
  epsilon = zeros(ng, ng); 
  for ind_nv = nv-nv_oper+1:nv
    if aqsFlag
      Mgvc = GWinfo.aqs{ind_nv}(:, nv+1:nv+nc_oper);
    else
      Mgvc = mtxel_sigma(ind_nv, GWinfo, options.Groundstate, ...
            (nv+1:nv+nc_oper));
    end
  	Mgvc = conj(Mgvc);
    
    Eden = ev(ind_nv) - ev(nv+1:nv+nc_oper);
    edenDRtmp = (-1.0 ./ (omega - Eden - mi*eta) ...
    + 1.0 ./ (omega + Eden + mi*eta));

    epsilon = epsilon + 2*Mgvc*(edenDRtmp.*Mgvc') / vol;
  end % for ind_nv

  if 1
    epsilon = eye(ng) - Dcoul * epsilon;
    epsilon = inv(epsilon);
    W = (eye(ng) - epsilon)*Dcoul / vol;
    if testflag1
      W_old = (I_eps_array(:, :, ifreq))*Dcoul/vol;
      out = norm(W - W_old);
      fprintf('ifreq = %d, epsilon difference = %.3e.\n', ifreq, out)
    end
    clear epsilon;
  end


  % Calculate M_{ni}' * W(\omega) *M_{ni} for both i and n are occupied states
  % Save it into MWM_occ
  for ibandener_count = 1:nv_ener
    ibandener = bandtocal_occ(ibandener_count);
    if 1
      Mgvc = mtxel_sigma(ibandener, GWinfo, options.Groundstate, ...
                        (nv-nv_oper+1):nv);
      Mgvc = conj(Mgvc);
      temp = W * Mgvc;
      for ibandoper_Mg = 1:nv_oper
        MWM_occ(ibandener_count, ibandoper_Mg, ifreq) ...
        = Mgvc(:, ibandoper_Mg)' * temp(:, ibandoper_Mg);
      end
    end
  end
  % Calculate for both i and n are unoccupied states.
  % Sae it into MWM_unocc
  for ibandener_count = 1:nc_ener 
    ibandener = bandtocal_unocc(ibandener_count);
    if 1
      Mgvc = mtxel_sigma(ibandener, GWinfo, options.Groundstate, ...
                        (nv+1):nv+nc_oper);
      Mgvc = conj(Mgvc);
      temp = W * Mgvc;
      for ibandoper_Mg = 1:nc_oper
        MWM_unocc(ibandener_count, ibandoper_Mg, ifreq) ...
        = Mgvc(:, ibandoper_Mg)' * temp(:, ibandoper_Mg);
      end
    end
  end
end% for ifreq

% if testflag1
%   for ifreq = 1:nfreq_real
%     norm(sres_elements(1:4, 1:4, ifreq) - MWM_occ(:, :, ifreq))
%     norm(sres_elements(5:8, 5:8, ifreq) - MWM_unocc(:, :, ifreq))
%   end
%   % sres_ele_old = sres_elements(ibandener_count, ibandoper_Mg, ifreq);
%   % sres_ele = MWM_occ(ibandener_count, ibandoper_Mg, ifreq)
%   % if (abs(sres_ele_old - sres(ibandener_count)) >= TOL_ZERO)
%   %   fprintf("n = %d, i = %d, sres not equal to sres_elements!", ibandener, ibandoper);
%   % end
% end

% Use the elements in MWM to calculate sres
% First, both occupied states.
occ_sign = -1;
for ibandener_count = 1:nv_ener 
  ibandener = bandtocal_occ(ibandener_count);
  for ibandoper = nv-nv_oper+1:nv
    if testflag1 
      coeffcount = 0;
    end
    ibandoper_Mg = ibandoper - nv+nv_oper;
    if flagry2ev
      %x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
      x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL*ry2ev;
    else
      x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
      % x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL/ry2ev;
    end
    if x >= 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      if testflag1
        coeff_old = coeff_real(ibandener_count, ibandoper_Mg, ifreq);
        if abs(coeff - coeff_old) >= TOL_ZERO
          fprintf("n = %d, i = %d, ifreq = %d, coeff not equal to coeff_real!\n", ...
                   ibandener, ibandoper, ifreq);
          fprintf("coeff = %.3e, coeff_old = %.3e.\n", coeff, coeff_old);
        end
        elements_new = MWM_occ(ibandener_count, ibandoper_Mg, ifreq);
        elements_old = sres_elements(ibandener_count, ibandoper_Mg, ifreq);
        if abs(elements_new - elements_old) >= TOL_ZERO
          fprintf("n = %d, i = %d, elements not equal to sres_elements!", ...
                   ibandener, ibandoper);
        end
      end 

      if abs(coeff) <= TOL_ZERO
        continue;
      end
      sres(ibandener_count) = sres(ibandener_count) ...
      - coeff * occ_sign * MWM_occ(ibandener_count, ibandoper_Mg, ifreq); 
      if testflag1
        coeffcount = coeffcount+coeff;
      end

      if testflag1
        if (abs(coeffcount - 1) >= TOL_ZERO)
          fprintf("n = %d, i = %d, coeff added up = %d, not 1!\n", ...
                  ibandener, ibandoper, coeffcount);
        end
      end
    end
  end
end

% Then both unoccupied states
occ_sign = 1;
for ibandener_count = 1:nc_ener 
  ibandener = bandtocal_unocc(ibandener_count);
  for ibandoper = nv+1:nv+nc_oper
    ibandoper_Mg = ibandoper - nv;
    if flagry2ev
      x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL*ry2ev;
%      x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    else
%       x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL/ry2ev;
      x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    end
    if x < 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      sres(ibandener_count+nv_ener) = sres(ibandener_count+nv_ener) ...
      - coeff * occ_sign * MWM_unocc(ibandener_count, ibandoper_Mg, ifreq); 
    end 
  end
end
sres = sres;
end % function
