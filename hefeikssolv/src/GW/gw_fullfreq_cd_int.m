function sint = gw_fullfreq_cd_int(ksinfo, options)
% This file will be final version for the full-frequency GW calculation.
% Generally speaking, user should not directly call this function, instead
% call this indirectly by calling through gwCalculation->gw_fullfreq_cd.

% Description
%     Calculate the self-energies using contour integral idea within the GW method.
%     This code only computes the integral on the imaginary axis,  

% Parameters
%   Input:
%     ksinfo: class @ksinfo, contains ground-state information.
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

% testflag1 is for debugging
% flag2 is for method (numerically equal, suggest to use flag2 = true.)
testflag1 = true;
testflag1 = false;
flag2 = true;
% flag2 = false;


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
mi = sqrt(-1);
bandtocal = nv-nv_ener+1:nv+nc_ener;
n_ener = length(bandtocal);

% Initialize other values
ksinfo.Z     = ksinfo.Z * sqrt(vol);
Z     = ksinfo.Z;
ev    = ksinfo.ev;
F     = ksinfo.F;
Vxc   = ksinfo.Vxc;
Dcoul = spdiags(ksinfo.coulG(:,4), 0, ng, ng);
aqsFlag = ~isempty(ksinfo.aqs);
Dcoul(1, 1) = ksinfo.coulG0;
sqrtDcoul = sqrt(Dcoul); 

% 1. Generate frequency sequences in imaginary axis.
startFrequencyGen = tic;
[grid_real, ~, grid_imag, coeff_imag_func] = freqgen(ksinfo, options);
nfreq_real = length(grid_real);
nfreq_imag = length(grid_imag);
% eta = options.GWCal.dBrdning;
eta = 0;
% eta = 1e-4;
timeFrequencyGen = toc(startFrequencyGen);
fprintf('Time for Generating Frequencies = %.3s.\n', timeFrequencyGen);

% Energies use unit 'ev' in this code.
% Since both KSSOLV_dft and frequency generating part use Ry as unit
% Change unit first.
ev = ev * ry2ev;
grid_real = grid_real * ry2ev;
grid_imag = grid_imag * ry2ev;
Dcoul = Dcoul * ry2ev;


sint = zeros(n_ener, 1);
if testflag1
  load intcoeff.mat
  load I_eps_array.mat
  load sint_n_i_omega.mat
  load aqslist.mat
  aqscount = 0;
end
for ifreq = 1:nfreq_imag
  omega = grid_imag(ifreq);
  coeff_func = coeff_imag_func{ifreq};
  % 2. Calculate W(\omega) for each frequency 
  epsilon = zeros(ng, ng);
  for ind_nv = nv-nv_oper+1:nv
    if aqsFlag
      Mgvc = ksinfo.aqs{ind_nv}(:, nv+1:nv+nc_oper);
    else
      Mgvc = mtxel_sigma(ind_nv, ksinfo, options.Groundstate, ...
            (nv+1:nv+nc_oper));
    end
  	Mgvc = conj(Mgvc);
    
    Eden = ev(ind_nv) - ev(nv+1:nv+nc_oper);
    edenDRtmp = (-1.0 ./ (omega - Eden - mi*eta) ...
    + 1.0 ./ (omega + Eden + mi*eta));
    
    epsilon = epsilon + 2*Mgvc*(edenDRtmp.*Mgvc') / vol;
  end % for ind_nv
  if ~flag2
    epsilon = eye(ng) - Dcoul * epsilon;
    epsilon = inv(epsilon);
    W = (eye(ng) - epsilon)*Dcoul / vol;
    if testflag1
      W_old = (I_eps_array(:, :, ifreq+nfreq_real))*Dcoul/vol;
      out = norm(W - W_old);
      fprintf('ifreq = %d, epsilon difference = %.3e.\n', ifreq, out)
      clear epsilon;
    end
  end
  if flag2
    epsilon = (Dcoul/vol)^(-1) - epsilon*vol;
    % Here epsilon is supposed to be an Hermite matrix.
    if testflag1
      epsilon_old = epsilon;
    end
    epsilon = tril(epsilon, -1) + tril(epsilon, -1)' + diag(real(diag(epsilon)));
    if testflag1
      norm(epsilon_old - epsilon)
    end
    [L, D] = ldl(epsilon); rsqrtD = diag(sqrt(diag(D)).^(-1));
    W = Dcoul/vol - (L')^(-1) * rsqrtD * rsqrtD * (L)^(-1);
    if testflag1
      W_old = (I_eps_array(:, :, ifreq+nfreq_real))*Dcoul/vol;
      out = norm(W - W_old);
      fprintf('ifreq = %d, epsilon difference = %.3e.\n', ifreq, out)
    end
  end
  % 3. For each band i to construct GW operators
  for ibandener_count = 1:length(bandtocal) % Iteration over n
    ibandener = bandtocal(ibandener_count);
    % Here, we calculate \sum_i coeff_n(ev_n-ev_i) W(\omega)\rho_ni^*\rho_ni.
    % Generate rho_nI
    if flag2
      Mgvc = mtxel_sigma(ibandener, ksinfo, options.Groundstate, ...
                        (nv-nv_oper+1):nv+nc_oper);
      Mgvc = conj(Mgvc);
      out_list = zeros(n_oper, 1);
      out_list = sum(Dcoul/vol*abs(Mgvc).^2)';
      Mgvc = rsqrtD * (L^(-1)*Mgvc); 
      out_list = out_list - sum(abs(Mgvc).^2)';
      clear Mgvc;
      for ibandoper = nv-nv_oper+1:nv+nc_oper  
        % x = (ev(ibandener) - ev(ibandoper))*ry2ev + TOL_SMALL;
        x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
        if abs(x) < TOL_SMALL
          coeff = coeff_func(TOL_SMALL);
        else
          coeff = coeff_func(x);
        end
        ibandoper_Mg = ibandoper - (nv-nv_oper);   
        if testflag1
          coeff_old = coeff_matrix(ibandener, ibandoper, ifreq);
          if abs(coeff_old - coeff) >= 1e-6
            fprintf("n = %d, i = %d, ifreq = %d,\n",...
            ibandener, ibandoper, ifreq);
            fprintf("coeff_old = %.6e, coeff = %.6e, diff = %.3e.\n", ...
            coeff_old, coeff, coeff_old - coeff);
          end
        end
        if (ibandener == 2 & ibandoper == 2)
          1;
        end
        % out = Mgvc(:, ibandoper_Mg)'*Mgvc(:, ibandoper_Mg); 
        out = out_list(ibandoper_Mg);
        if testflag1
          out_old = sint_out(ibandener, ibandoper, ifreq);
          fprintf("n = %d, i = %d, ifreq = %d,\n",...
          ibandener, ibandoper, ifreq);
          fprintf("func_int = %.6e, func = %.6e, diff = %.3e.\n", ...
          out_old, out, out_old - out);
        end
        sint(ibandener_count) = sint(ibandener_count) + coeff*out;
      end
    end % if 0
    if ~flag2
      Mgvc = mtxel_sigma(ibandener, ksinfo, options.Groundstate, ...
                        (nv-nv_oper+1):nv+nc_oper);
      Mgvc = conj(Mgvc);
      for ibandoper = nv-nv_oper+1:nv+nc_oper  
        % x = (ev(ibandener) - ev(ibandoper))*ry2ev + TOL_SMALL;
        x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
        if testflag1
          coeff_old = coeff_matrix(ibandener, ibandoper, ifreq);
        end
        if abs(x) < TOL_SMALL
          coeff = coeff_func(TOL_SMALL);
        else
          coeff = coeff_func(x);
        end
        ibandoper_Mg = ibandoper - (nv-nv_oper);   
        if testflag1
          if abs(coeff_old - coeff) >= 1e-6
            fprintf("n = %d, i = %d, ifreq = %d,\n",...
            ibandener, ibandoper, ifreq);
            fprintf("coeff_old = %.6e, coeff = %.6e, diff = %.3e.\n", ...
            coeff_old, coeff, coeff_old - coeff);
          end
          if ifreq == 1
            aqscount = aqscount+1;
            if (norm(aqslist(:, aqscount) - Mgvc(:, ibandoper_Mg)) >= 1e-8)
              fprintf("n = %d, i = %d, ifreq = %d,\n",...
              ibandener, ibandoper, ifreq);
              fprintf("aqstemp is different!");
            end
          end
        end
        if (ibandener == 2 & ibandoper == 2)
          1;
        end
        out = Mgvc(:, ibandoper_Mg)'*W*Mgvc(:, ibandoper_Mg); 
        if testflag1
          out_old = sint_out(ibandener, ibandoper, ifreq);
          fprintf("n = %d, i = %d, ifreq = %d,\n",...
          ibandener, ibandoper, ifreq);
          fprintf("func_int = %.6e, func = %.6e, diff = %.3e.\n", ...
          out_old, out, out_old - out);
        end
        sint(ibandener) = sint(ibandener) + coeff*out;
      end
    end % if 0
  end
end 

% sint = sint / pi * ry2ev;
sint = sint / pi;
[real(sint), imag(sint)];

timeintegral = toc(startintgral);
fprintf("Time to integral on Imaginary axis = %.2d.\n", timeintegral);
end