function sint = gw_fullfreq_cd_int_ISDF(GWinfo, options)
% This file will be final version for the full-frequency GW calculation.
% Generally speaking, user should not directly call this function, instead
% call this indirectly by calling through gwCalculation->gw_fullfreq_cd.

% Description
%     Calculate the self-energies using contour integral idea within the GW method.
%     This code only computes the integral on imaginary axis.
%     This code applies ISDF algorithm. 

% Parameters
%   Input:
%     GWinfo: class @GWinfo, contains ground-state information.
%     options: class @GWOptions, contains necessary parameters for the calculation.
%   Output:
%     sint: the integral on imaginary axis, size nbandener*1.

% Structure 
%   The code is organized as followed:
%     0. Initialize inputs.
%     1. Do ISDF if no ISDFoutput in current folder, else read the files.
%     2. Generate the qudrature weights and the quadrature node on the imaginary axis.
%        For each frequencies omega
%     3    Calculate K(omega), where K is a frequency-dependent kernel,
%          and do the LDLT decomposition. 
%          For each band to calculate n
%     4.     Calculate <nn'|W(omega)|nn'> using K(omega) for all n'. 
%     5.     Multiply the coefficient and add the result to sint  
%          end 
%        end 

% testflag1 is for debugging
% flag2 is for method (numerically equal, suggest to use flag2 = true.)
flag2 = 1;


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
Z     = GWinfo.Z;
ev    = GWinfo.ev;
gvec = GWinfo.gvec;
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
Dcoul(1, 1) = GWinfo.coulG0;

% Energies use unit 'ev' in this code.
% Since both KSSOLV_dft and frequency generating part use Ry as unit
% Change unit first.
ev = ev * ry2ev;
Dcoul = Dcoul * ry2ev;
sqrtDcoul = sqrt(Dcoul); 


% 1. Generate the qudrature weights and the quadrature node
%    on the imaginary axis.
startFrequencyGen = tic;
[~, ~, grid_imag, coeff_imag_func] = freqgen(GWinfo, options);
nfreq_imag = length(grid_imag);
% eta = options.GWCal.dBrdning;
eta = 0;
% eta = 1e-4;
timeFrequencyGen = toc(startFrequencyGen);
fprintf('Time for Generating Frequencies = %.3s.\n', timeFrequencyGen);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start ISDF
% Generally, we have relations 
%     Mr{xx} = xxrzeta_mu * Mrxx_mu;
startISDF = tic;

vcrank_mu = ceil(sqrt(nv_oper*nc_oper) * options.ISDFCauchy.vcrank_ratio);
vsrank_mu = ceil(sqrt(nv_oper*n_ener)  * options.ISDFCauchy.vsrank_ratio);
ssrank_mu = ceil(sqrt(n_ener *n_ener)  * options.ISDFCauchy.ssrank_ratio);



psir = zeros(nr, nv+nc);
for iband = 1:nv+nc
  fftbox1 = put_into_fftbox(Z(:, iband), gvec.idxnz, gvec.fftgrid);
  fftbox1 = gvec.nfftgridpts / GWinfo.vol * do_FFT(fftbox1, gvec.fftgrid, 1);
  psir(:, iband) = reshape(fftbox1, gvec.nfftgridpts, []);
end
vcrzeta_mu = zeros(nr, vcrank_mu);
vcgzeta_mu = zeros(ng, vcrank_mu);
optionsISDF = options.ISDFCauchy;


startVC = tic;
optionsISDF.isdfoptions.rank = vcrank_mu;
% if ~flag3
%   [vcrzeta_mu, vcind_mu] = isdf(conj(psir(:,nv-nv_oper+1:nv)), ...
%                                (psir(:,nv+1:nv+nc_oper)), optionsISDF);
%   for i = 1:vcrank_mu
%     fftbox1 = reshape(vcrzeta_mu(:, i), gvec.fftgrid);
%     fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * GWinfo.vol;
%     vcgzeta_mu(:, i) = get_from_fftbox(gvec.idxnz, fftbox1, gvec.fftgrid);
%   end
%   clear vcrzeta_mu;
% else
%   vcind_mu = isdf_indices(conj(psir(:,nv-nv_oper+1:nv)), ...
%                                (psir(:,nv+1:nv+nc_oper)), optionsISDF);
%   vcgzeta_mu = isdf_kernelg(conj(psir(:,nv-nv_oper+1:nv)), ...
%               (psir(:,nv+1:nv+nc_oper)), vcind_mu, gvec, vol);
% end
[vcind_mu, vcgzeta_mu] = isdf_main('vc', psir, nv-nv_oper+1:nv, ...
        nv+1:nv+nc_oper, gvec, vol, optionsISDF);
vcgzeta_mu = conj(vcgzeta_mu);
timeforVC = toc(startVC);


startVS = tic;
vsrzeta_mu = zeros(nr, vsrank_mu);
vsgzeta_mu = zeros(ng, vsrank_mu);
optionsISDF.isdfoptions.rank = vsrank_mu;
[vsind_mu, vsgzeta_mu] = isdf_main('vs', psir, nv-nv_oper+1:nv, ...
        nv-nv_ener+1:nv+nc_ener, gvec, vol, optionsISDF);
vsgzeta_mu = conj(vsgzeta_mu);
timeforVS = toc(startVS);


startSS = tic;
ssrzeta_mu = zeros(nr, vsrank_mu);
ssgzeta_mu = zeros(ng, ssrank_mu);
optionsISDF.isdfoptions.rank = ssrank_mu;
[ssind_mu, ssgzeta_mu] = isdf_main('ss', psir, nv-nv_oper+1:nv+nc_oper, ...
        nv-nv_ener+1:nv+nc_ener, gvec, vol, optionsISDF);
ssgzeta_mu = conj(ssgzeta_mu);
timeforSS = toc(startSS);
timeforISDF = toc(startISDF);

% ssKernel = ssgzeta_mu' * Dcoul * ssgzeta_mu;
vcVvc = vcgzeta_mu' * Dcoul * vcgzeta_mu / vol;
vcVnn = vcgzeta_mu' * Dcoul * ssgzeta_mu / vol;
% cholvcVvn = chol(vcVnn, 'lower');

sint = zeros(n_ener, 1);
for ifreq = 1:nfreq_imag
  % 3. Calculate K(omega) for each frequency, do the LDLT decomposition 
  if ~flag2 
    omega = grid_imag(ifreq);
    coeff_func = coeff_imag_func{ifreq};
    epsilon = zeros(ng, ng);
    for ind_nv = nv-nv_oper+1:nv
      Mgvc = (psir(vcind_mu, ind_nv)) .* conj(psir(vcind_mu, nv+1:nv+nc_oper)); 
      Eden = ev(ind_nv) - ev(nv+1:nv+nc_oper);
      edenDRtmp = (-1.0 ./ (omega - Eden - mi*eta) ...
      + 1.0 ./ (omega + Eden + mi*eta));
      epsilon = epsilon + 2*vcgzeta_mu * Mgvc*diag(edenDRtmp)*Mgvc'*vcgzeta_mu' / vol;
    end % for ind_nv
    epsilon = eye(ng) - Dcoul * epsilon;
    epsilon = inv(epsilon);
    W = (eye(ng) - epsilon)*Dcoul / vol;
    ssWss = ssgzeta_mu' * W * ssgzeta_mu;
  end
  if flag2
    omega = grid_imag(ifreq);
    coeff_func = coeff_imag_func{ifreq};
    epsKernel = zeros(vcrank_mu, vcrank_mu);
    for ind_nv = nv-nv_oper+1:nv
      Mgvc = (psir(vcind_mu, ind_nv)) .* conj(psir(vcind_mu, nv+1:nv+nc_oper)); 
      Eden = ev(ind_nv) - ev(nv+1:nv+nc_oper);
      edenDRtmp = (-1.0 ./ (omega - Eden - mi*eta) ...
      + 1.0 ./ (omega + Eden + mi*eta));
      epsKernel = epsKernel + Mgvc*diag(edenDRtmp)*Mgvc';
    end % for ind_nv 
    epsKernel = inv(epsKernel)/2 - vcVvc; 
    % Here epsKernel is supposed to be an Hermite matrix.
    epsKernel = tril(epsKernel, -1) + tril(epsKernel, -1)' + diag(real(diag(epsKernel)));
    [LKernel, DKernel] = ldl(epsKernel);
    dKernel = diag(DKernel);
    % rsqrtD = diag(sqrt(diag(D)).^(-1));
  end
  % Now, W is supposed to be VP_vc (LDL^T)^(-1) P_vcV
  % What we care about is
  % \sum_{i, mu1, mu2, omega'} coeff(n, i; omega') ...
  %     * (psi_i.*psi_n)' * P * W(omega') * P * (psi_i.*psi_n), 
  % and P_nn'*W*P_nn = - (L\vcVnn)' * D^(-1) * (L\vcVvn)
  
  % So we construct P_nn'*W(omega')*P_nn first.
  % Since what we care is a quadric form, we only need "half of" it
  if flag2
    right_nnWnn = LKernel \ vcVnn;
  end

  % 4. For each band, calculate <nn'|W(omega)|nn'> using K(\omega)
  for ibandener_count = 1:length(bandtocal) % Iteration over n
    ibandener = bandtocal(ibandener_count);
    % Here, we calculate \sum_i coeff_n(ev_n-ev_i) W(\omega)\rho_ni^*\rho_ni.
    % Generate rho_nI
    if flag2
      Mgvc = psir(ssind_mu, ibandener) ...
               .* conj(psir(ssind_mu, nv-nv_oper+1:nv+nc_oper)); 
      Mgvc = right_nnWnn*Mgvc; 
      out_list = - sum(dKernel.^(-1) .* abs(Mgvc).^2);
      clear Mgvc;
      % Mgvc = mi* + sqrt(Dcoul/vol)\Mgvc;
      for ibandoper = nv-nv_oper+1:nv+nc_oper  
        x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
        if abs(x) < TOL_SMALL
          coeff = coeff_func(TOL_SMALL);
        else
          coeff = coeff_func(x);
        end
        ibandoper_Mg = ibandoper - (nv-nv_oper);   
        % out = Mgvc(:, ibandoper_Mg)'*Mgvc(:, ibandoper_Mg); 
        out = out_list(ibandoper_Mg);
        sint(ibandener_count) = sint(ibandener_count) + coeff*out;
      end
    end % if 0
    if ~flag2
      C_in_n = conj(psir(ssind_mu, ibandener))... 
               .* psir(ssind_mu, nv-nv_oper+1:nv+nc_oper); 
      C_in_n = conj(C_in_n);
      for ibandoper = nv-nv_oper+1:nv+nc_oper  
        x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
        if abs(x) < TOL_SMALL
          coeff = coeff_func(TOL_SMALL);
        else
          coeff = coeff_func(x);
        end
        ibandoper_Mg = ibandoper - (nv-nv_oper);   
        out = C_in_n(:, ibandoper_Mg)'*ssWss*C_in_n(:, ibandoper_Mg); 
        sint(ibandener) = sint(ibandener) + coeff*out;
      end
    end % if 0
  end
end 

sint = sint / pi;
% [real(sint), imag(sint)]

timeintegral = toc(startintgral);
fprintf("Time to integral on Imaginary axis = %.2d.\n", timeintegral);
end
