function sint = gw_fullfreq_cd_int(GWinfo, options)
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

% flag2 is for method (numerically equal, suggest to use flag2 = true.)
flag2 = false;
testflag1 = false;
testflag2 = false;
% testflag2 = true;
ii = sqrt(-1);
ind_bgw2ks = options.Groundstate.ind_bgw2ks;

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
GWinfo.Z     = GWinfo.Z * sqrt(vol);
ev    = GWinfo.ev;
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
aqsFlag = ~isempty(GWinfo.aqs);
Dcoul(1, 1) = GWinfo.coulG0;

% 1. Generate frequency sequences in imaginary axis.
startFrequencyGen = tic;
[grid_real, ~, grid_imag, coeff_imag_func] = freqgen(GWinfo, options);
nfreq_real = length(grid_real);
nfreq_imag = length(grid_imag);
eta = 0;
% eta = 1e-4;
timeFrequencyGen = toc(startFrequencyGen);
fprintf('Time for Generating Frequencies = %.3s.\n', timeFrequencyGen);

% Energies use unit 'ev' in this code.
% Since both KSSOLV_dft and frequency generating part use Ry as unit
% Change unit first.
ev = ev * ry2ev;
Dcoul = Dcoul * ry2ev;

% ! testperpose

sint = zeros(n_ener, 1);

if testflag1
  W_V = readmatrix('W_V');
  W_V = W_V(:, 1) + ii*W_V(:, 2);
  W_V = reshape(W_V, ng, ng, []);
  for i = 1:nfreq_real + nfreq_imag
    W_V(:, :, i) = W_V(ind_bgw2ks, ind_bgw2ks, i);
  end
end

for ifreq = 1:nfreq_imag
  omega = grid_imag(ifreq);
  coeff_func = coeff_imag_func{ifreq};
  % 2. Calculate W(\omega) for each frequency 
  epsilon = zeros(ng, ng);
  Dcoul(1,1) = 0.0;
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
    
%    epsilon = epsilon + 2*Mgvc*(edenDRtmp.*Mgvc');
    epsilon = epsilon + 2*Mgvc*(edenDRtmp.*Mgvc') / vol;
  end % for ind_nv
  chinorm = norm(epsilon, 'fro') * ry2ev;
  fprintf('ifreq = %6d, chinorm = %12.6f\n', ifreq+nfreq_real, chinorm);
  if ~flag2
    epsilon = eye(ng) - Dcoul * epsilon;
    chinorm = norm(epsilon, 'fro');
    fprintf('ifreq = %6d, epsnorm = %12.6f\n', ifreq+nfreq_real, chinorm);
    epsilon = inv(epsilon);
    Dcoul(1, 1) = GWinfo.coulG0 * ry2ev;
    W = (eye(ng) - epsilon) * Dcoul / vol;
  end
  if flag2
    epsilon = (Dcoul/vol)^(-1) - epsilon*vol;
    % Here epsilon is supposed to be an Hermite matrix.
    epsilon = tril(epsilon, -1) + tril(epsilon, -1)' + diag(real(diag(epsilon)));
    [L, D] = ldl(epsilon); rsqrtD = diag(sqrt(diag(D)).^(-1));
    W = Dcoul/vol - (L')^(-1) * rsqrtD * rsqrtD * (L)^(-1);
  end
  % 3. For each band i to construct GW operators
  for ibandener_count = 1:length(bandtocal) % Iteration over n
    ibandener = bandtocal(ibandener_count);
    % Here, we calculate \sum_i coeff_n(ev_n-ev_i) W(\omega)\rho_ni^*\rho_ni.
    % Generate rho_nI
    if flag2
      Mgvc = mtxel_sigma(ibandener, GWinfo, options.Groundstate, ...
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
        if (ibandener == 2 & ibandoper == 2)
          1;
        end
        % out = Mgvc(:, ibandoper_Mg)'*Mgvc(:, ibandoper_Mg); 
        out = out_list(ibandoper_Mg);
        sint(ibandener_count) = sint(ibandener_count) + coeff*out;
      end
    end % if 0
    if ~flag2
      Mgvc = mtxel_sigma(ibandener, GWinfo, options.Groundstate, ...
                        (nv-nv_oper+1):nv+nc_oper);
      Mgvc = conj(Mgvc);
      if testflag2
        aqs_file = ['aqs_', num2str(ibandener_count)];
        aqs = readmatrix(aqs_file);
        aqs = aqs(:, 1) + ii*aqs(:, 2);
        aqs = reshape(aqs, ng, []);
        aqs(:, :) = aqs(ind_bgw2ks, :);
        fprintf("Fnorm(aqs(%d)) = %12.6f\n", ibandener_count, norm(aqs, 'fro'));
        fprintf("Fnorm(Mg(%d)) = %12.6f\n", ibandener_count, norm(aqs, 'fro'));
      end
      for ibandoper = nv-nv_oper+1:nv+nc_oper  
        % x = (ev(ibandener) - ev(ibandoper))*ry2ev + TOL_SMALL;
        x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
        if abs(x) < TOL_SMALL
          coeff = coeff_func(TOL_SMALL);
        else
          coeff = coeff_func(x);
        end
%         fprintf("iener = %6d, ioper = %6d, coeff = %12.6f\n", ibandoper, ibandener, coeff); 
        ibandoper_Mg = ibandoper - (nv-nv_oper);   
        if (ibandener == 2 & ibandoper == 2)
          1;
        end
        out = Mgvc(:, ibandoper_Mg)'*W*Mgvc(:, ibandoper_Mg); 
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
