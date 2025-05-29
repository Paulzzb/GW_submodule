function sres = gw_fullfreq_cd_res_ISDF(GWinfo, options)
% This file will be final version for the full-frequency GW calculation.
% Generally speaking, user should not directly call this function, instead
% call this indirectly by calling through gwCalculation->gw_fullfreq_cd.

% Description
%     Calculate the self-energies using contour integral idea within the GW method.
%     This code only computes the summation of the residual.
%     This code applies ISDF algorithm. 

% Parameters
%   Input:
%     GWinfo: class @GWinfo, contains ground-state information.
%     options: class @GWOptions, contains necessary parameters for the calculation.
%   Output:
%     sres: the summation of the residual, size nbandener*1.

% Structure 
%   The code is organized as followed:
%     0. Initialize inputs.
%     1. Do ISDF if no ISDFoutput in current folder, else read the files.
%     2. Generate the qudrature weights and the quadrature node on the real axis.
%        For each frequencies omega
%     3. Generate the qudrature weights and the quadrature node on the imaginary axis.
%        For each frequencies omega
%     4    Calculate K(omega), where K is a frequency-dependent kernel,
%          and do the LDLT decomposition. 
%          For each band to calculate n
%     5.     Calculate <nn'|W(omega)|nn'> using K(omega) for all n'. 
%     6.     Multiply the coefficient and add the result to sint  
%          end 
%        end 


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
Z        = GWinfo.Z * sqrt(vol);
gvec = GWinfo.gvec;
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
Dcoul(1, 1) = GWinfo.coulG0;

% We use unit as ev, while usual input is Ry.
ev    = GWinfo.ev;
Dcoul = Dcoul * ry2ev;


% 1. Generate frequency sequences in real axis and imaginary axis.
startFrequencyGen = tic;
[grid_real, coeff_real_func, ~, ~] = freqgen(GWinfo, options);
nfreq_real = length(grid_real);
% eta = options.GWCal.dBrdning / ry2ev;
eta = 0;
timeFrequencyGen = toc(startFrequencyGen);
fprintf('Time for Generating Frequencies = %.3s.\n', timeFrequencyGen);


% ISDF here
vcrank_mu = ceil(sqrt(nv_oper*nc_oper) * options.ISDFCauchy.vcrank_ratio);
vsrank_mu = ceil(sqrt(nv_oper*n_ener)  * options.ISDFCauchy.vsrank_ratio);
ssrank_mu = ceil(sqrt(n_ener *n_ener)  * options.ISDFCauchy.ssrank_ratio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start ISDF
% Generally, we have relations 
%     Mr{xx} = xxrzeta_mu * Mrxx_mu;
startISDF = tic;

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





sres = zeros(n_ener, 1);

MWM_occ = zeros(nv_ener, nv_oper, nfreq_real);
MWM_unocc = zeros(nc_ener, nc_oper, nfreq_real);
vcVvc = vcgzeta_mu' * Dcoul * vcgzeta_mu / vol;
vcVnn = vcgzeta_mu' * Dcoul * ssgzeta_mu / vol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Calculation


for ifreq = 1:nfreq_real
  omega = grid_real(ifreq);
  epsilon = zeros(vcrank_mu, vcrank_mu); 
  for ind_nv = nv-nv_oper+1:nv
    Mgvc = psir(vcind_mu, ind_nv) .* conj(psir(vcind_mu, nv+1:nv+nc_oper)); 
    Eden = ev(ind_nv) - ev(nv+1:nv+nc_oper);
    edenDRtmp = (-1.0 ./ (omega - Eden - mi*eta) ...
    + 1.0 ./ (omega + Eden + mi*eta));
    epsilon = epsilon + Mgvc*diag(edenDRtmp)*Mgvc';
  end % for ind_nv
  epsKernel = inv(epsilon)/2 - vcVvc; 
  clear epsilon
  WKernel = -inv(epsKernel);
  clear epsKernel;

  % Now, W is supposed to be V*P_vc * epsKernel * P_vc'*V

  % Calculate M_{ni}' * W(\omega) *M_{ni} for both i and n are occupied states
  % Save it into MWM_occ
  for ibandener_count = 1:nv_ener
    ibandener = bandtocal_occ(ibandener_count);
    Mg = psir(ssind_mu, ibandener) .* conj(psir(ssind_mu, (nv-nv_oper+1):nv));
    Mg = vcVnn * Mg;
    temp = WKernel*Mg;
    for ibandoper_Mg = 1:nv_oper
      if (MWM_occ(ibandener_count, ibandoper_Mg, ifreq) > 0)
        MWM_occ(ibandener_count, ibandoper_Mg, ifreq) ...
        = Mg(:, ibandoper_Mg)' * temp(:, ibandoper_Mg);
      end
    end
  end
  % Calculate for both i and n are unoccupied states.
  % Save it into MWM_unocc
  for ibandener_count = 1:nc_ener 
    ibandener = bandtocal_unocc(ibandener_count);
    Mg = psir(ssind_mu, ibandener) .* conj(psir(ssind_mu, nv+1:nv+nc_oper));
    Mg = vcVnn * Mg;
    temp = WKernel * Mg;
    for ibandoper_Mg = 1:nc_oper
      if (MWM_unocc(ibandener_count, ibandoper_Mg, ifreq) > 0)
        MWM_unocc(ibandener_count, ibandoper_Mg, ifreq) ...
        = Mg(:, ibandoper_Mg)' * temp(:, ibandoper_Mg);
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
      sres(ibandener_count) = sres(ibandener_count) ...
      - coeff * occ_sign * MWM_occ(ibandener_count, ibandoper_Mg, ifreq); 
    end
  end
end

% Then both unoccupied states
occ_sign = 1;
for ibandener_count = 1:nc_ener 
  ibandener = bandtocal_unocc(ibandener_count);
  for ibandoper = nv+1:nv+nc_oper
    ibandoper_Mg = ibandoper - nv;
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
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
