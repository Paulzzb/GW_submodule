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
%     2. Calculate pattern, i.e., for which set of (n, n') we need to 
%        calculate <nn'|W(omega)|nn'>.
%     3    Calculate W(omega)
%     4.   Calculate <nn'|W(omega)|nn'> for all n and n' both occupied
%          and both unoccupied, then save it. 
%          For each n and n' both occupied and unoccupied
%     5.     Check if it is residual, add to sres.
%          end 
%        end 


testflag1 = false;
testflag2 = false;


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
% We use unit as ev, while usual input is Ry.
GWinfo.Z = GWinfo.Z * sqrt(vol);
Z        = GWinfo.Z * sqrt(vol);
ev       = GWinfo.ev * ry2ev;
Dcoul    = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
aqsFlag  = ~isempty(GWinfo.aqs);
Dcoul(1, 1) = GWinfo.coulG0;
Dcoul    = Dcoul * ry2ev;


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


sres = zeros(n_ener, 1);
% if testflag1
%   W_V = readmatrix('W_V');
%   W_V = W_V(:, 1) + ii*W_V(:, 2);
%   W_V = reshape(W_V, ng, ng, []);
%   for i = 1:nfreq_real + nfreq_imag
%     W_V(:, :, i) = W_V(ind_bgw2ks, ind_bgw2ks, i);
%   end
% end




% Place to save the <nn'|W(omega)|nn'>, MWM_occ for two occupied,
% and MWM_unocc for two unoccupied.
MWM_occ = zeros(nv_ener, nv_oper, nfreq_real);
MWM_unocc = zeros(nc_ener, nc_oper, nfreq_real);

% Calculate pattern first, deciding which part to calculate
% After that, MWM_occ/MWM_unocc contains 1 or 0.
% We only need values for place 1.
for ibandener_count = 1:nv_ener 
  ibandener = bandtocal_occ(ibandener_count);
  for ibandoper = nv-nv_oper+1:nv
    ibandoper_Mg = ibandoper - nv+nv_oper;
    % x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL*ry2ev;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    if x >= 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      if abs(coeff) <= TOL_ZERO
        continue;
      end
      MWM_occ(ibandener_count, ibandoper_Mg, ifreq) = 1;
    end
  end
end

% Then both unoccupied states
for ibandener_count = 1:nc_ener 
  ibandener = bandtocal_unocc(ibandener_count);
  for ibandoper = nv+1:nv+nc_oper
    ibandoper_Mg = ibandoper - nv;
    % x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL*ry2ev;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    if x < 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      if abs(coeff) > TOL_ZERO
        MWM_unocc(ibandener_count, ibandoper_Mg, ifreq) = 1;
      end
    end 
  end
end



for ifreq = 1:nfreq_real
  omega = grid_real(ifreq);
  epsilon = zeros(ng, ng); 
  % Dcoul(1, 1) = 0;
  for ind_nv = nv-nv_oper+1:nv
    if aqsFlag
      Mgvc = GWinfo.aqs{ind_nv}(:, nv+1:nv+nc_oper);
    else
      Mgvc = mtxel_sigma(ind_nv, GWinfo, options.Groundstate, ...
            (nv+1:nv+nc_oper));
    end
    Mgvc = conj(Mgvc);

    
    Eden = ev(ind_nv) - ev(nv+1:nv+nc_oper);
    % edenDRtmp = (-1.0 ./ (omega - Eden - mi*eta) ...
    % + 1.0 ./ (omega + Eden + mi*eta));
    edenDRtmp = (1.0 ./ (Eden - omega) ...
    + 1.0 ./ (Eden + omega));

    epsilon = epsilon + 2*Mgvc*(edenDRtmp.*Mgvc') / vol;
  end % for ind_nv
  fprintf('ifreq = %3d, norm(chi) = %12.6f\n', ifreq, norm(epsilon, 'fro')*ry2ev);

  if 1
    epsilon = eye(ng) - Dcoul * epsilon;
    epsilon = inv(epsilon);
    fprintf('ifreq = %3d, norm(epsilon) = %12.6f\n', ifreq, norm(epsilon, 'fro'));
    Dcoul(1, 1) = GWinfo.coulG0*ry2ev;
    W = (eye(ng) - epsilon)*Dcoul / vol;
    if testflag1
      fprintf("Fnorm(W_V(%d)) = %12.6f\n", ibandener_count, norm(W_V(:, :, ifreq), 'fro'));
      fprintf("Fnorm(W_me(%d)) = %12.6f\n", ibandener_count, norm(W(:, :, ifreq), 'fro'));
    end
    clear epsilon;
  end


  % Calculate M_{ni}' * W(\omega) *M_{ni} for both i and n are occupied states
  % Save it into MWM_occ
  for ibandener_count = 1:nv_ener
    ibandener = bandtocal_occ(ibandener_count);
    Mgvc = mtxel_sigma(ibandener, GWinfo, options.Groundstate, ...
                      (nv-nv_oper+1):nv);
    Mgvc = conj(Mgvc);

    temp = W * Mgvc;
    for ibandoper_Mg = 1:nv_oper
      if (MWM_occ(ibandener_count, ibandoper_Mg, ifreq) > 0)
        MWM_occ(ibandener_count, ibandoper_Mg, ifreq) ...
        = Mgvc(:, ibandoper_Mg)' * temp(:, ibandoper_Mg);
      end
    end
  end
  % Calculate for both i and n are unoccupied states.
  % Save it into MWM_unocc
  for ibandener_count = 1:nc_ener 
    ibandener = bandtocal_unocc(ibandener_count);
    Mgvc = mtxel_sigma(ibandener, GWinfo, options.Groundstate, ...
                      (nv+1):nv+nc_oper);
    Mgvc = conj(Mgvc);
    temp = W * Mgvc;
    for ibandoper_Mg = 1:nc_oper
      if (MWM_unocc(ibandener_count, ibandoper_Mg, ifreq) > 0)
        MWM_unocc(ibandener_count, ibandoper_Mg, ifreq) ...
        = Mgvc(:, ibandoper_Mg)' * temp(:, ibandoper_Mg);
      end
    end
  end
end% for ifreq


% Use the elements in MWM to calculate sres
% First, both occupied states.
occ_sign = -1;
for ibandener_count = 1:nv_ener 
  ibandener = bandtocal_occ(ibandener_count);
  for ibandoper = nv-nv_oper+1:nv
    ibandoper_Mg = ibandoper - nv+nv_oper;
    % x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL*ry2ev;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    
    if x >= 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      if abs(coeff) > TOL_ZERO
        sres(ibandener_count) = sres(ibandener_count) ...
        - coeff * occ_sign * MWM_occ(ibandener_count, ibandoper_Mg, ifreq); 
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
    % x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL*ry2ev;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
    if x < 0
      continue;
    end
    x = abs(x);
    for ifreq = 1:nfreq_real
      coeff = coeff_real_func{ifreq}(x);
      if abs(coeff) > TOL_ZERO
        sres(ibandener_count+nv_ener) = sres(ibandener_count+nv_ener) ...
        - coeff * occ_sign * MWM_unocc(ibandener_count, ibandoper_Mg, ifreq); 
      end
    end 
  end
end
sres = sres;
end % function
