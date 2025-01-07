function sres = gw_fullfreq_cd_res(GWinfo, options)
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
%     2. Generate the qudrature weights and the quadrature node on the imaginary axis.
%        For each frequencies omega
%     3    Calculate K(omega), where K is a frequency-dependent kernel,
%          and do the LDLT decomposition. 
%          For each band to calculate n
%     4.     Calculate <nn'|W(omega)|nn'> using K(omega) for all n'. 
%     5.     Multiply the coefficient and add the result to sint  
%          end 
%        end 


testflag1 = true;
testflag1 = false;
testflag2 = true;
% testflag2 = false;
flagry2ev = true;
flagry2ev = false;


% Start the GW calculation
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
bandtocal_occ = find(bandtocal <= nv);
bandtocal_unocc = find(bandtocal > nv);
n_ener = length(bandtocal);
nv_ener = length(bandtocal_occ);
nc_ener = length(bandtocal_unocc);

% Initialize other values
% Energies use unit 'ev' in this code.
% Since both KSSOLV_dft and frequency generating part use Ry as unit
if flagry2ev
  GWinfo.ev = GWinfo.ev * ry2ev;
  GWinfo.Vxc = GWinfo.Vxc * ry2ev;
  options.GWCal.dFreqBrd = options.GWCal.dFreqBrd * ry2ev;
  options.GWCal.dFreqGrid = options.GWCal.dFreqGrid * ry2ev;
end
% 
GWinfo.Z     = GWinfo.Z * sqrt(vol);
Z     = GWinfo.Z;
ev    = GWinfo.ev;
Vxc   = GWinfo.Vxc;
gvec = GWinfo.gvec;
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
aqsFlag = ~isempty(GWinfo.aqs);
Dcoul(1, 1) = GWinfo.coulG0;
if flagry2ev
  Dcoul = Dcoul * ry2ev;
end


% 1. Generate frequency sequences in real axis and imaginary axis.
startFrequencyGen = tic;
[grid_real, coeff_real_func, grid_imag, coeff_imag_func] = freqgen(GWinfo, options);
nfreq_real = length(grid_real);
nfreq_imag = length(grid_imag);
% DOUBT!!!
if flagry2ev
  eta = options.GWCal.dBrdning;
else
  eta = options.GWCal.dBrdning / ry2ev;
  eta = 0;
end
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
% [vcrzeta_mu, vcind_mu] = isdf_gw(conj(psir(:,nv-nv_oper+1:nv)), ...
%                               (psir(:,nv+1:nv+nc_oper)), vcrank_mu);
% [vcrzeta_mu, vcind_mu] = isdf(conj(psir(:,nv-nv_oper+1:nv)), ...
%                              (psir(:,nv+1:nv+nc_oper)), optionsISDF);
% for i = 1:vcrank_mu
%   fftbox1 = reshape(vcrzeta_mu(:, i), gvec.fftgrid);
%   fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * GWinfo.vol;
%   vcgzeta_mu(:, i) = get_from_fftbox(gvec.idxnz, fftbox1, gvec.fftgrid);
% end
% vcgzeta_mu = conj(vcgzeta_mu);
% clear vcrzeta_mu;
% vcind_mu = isdf_indices(conj(psir(:,nv-nv_oper+1:nv)), ...
%                              (psir(:,nv+1:nv+nc_oper)), optionsISDF);
% vcgzeta_mu = isdf_kernelg(conj(psir(:,nv-nv_oper+1:nv)), ...
%             (psir(:,nv+1:nv+nc_oper)), vcind_mu, gvec, vol);
[vcind_mu, vcgzeta_mu] = isdf_main('vc', conj(psir(:,nv-nv_oper+1:nv)), ...
        (psir(:,nv+1:nv+nc_oper)), gvec, vol, optionsISDF);
vcgzeta_mu = conj(vcgzeta_mu);
timeforVC = toc(startVC);


startVS = tic;
vsrzeta_mu = zeros(nr, vsrank_mu);
vsgzeta_mu = zeros(ng, vsrank_mu);
optionsISDF.isdfoptions.rank = vsrank_mu;
% [vsrzeta_mu, vsind_mu] = isdf(conj(psir(:,nv-nv_oper+1:nv)), ...
                              % (psir(:,nv-nv_ener+1:nv+nc_ener)), optionsISDF);
% [vsrzeta_mu, vsind_mu] = isdf_gw(conj(psir(:,nv-nv_oper+1:nv)), ...
%                               (psir(:,nv-nv_ener+1:nv+nc_ener)), vsrank_mu);
% for i = 1:vsrank_mu
%   fftbox1 = reshape(vsrzeta_mu(:, i), gvec.fftgrid);
%   fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * GWinfo.vol;
%   vsgzeta_mu(:, i) = get_from_fftbox(gvec.idxnz, fftbox1, gvec.fftgrid);
% end
% vsgzeta_mu = conj(vsgzeta_mu);
% clear vsrzeta_mu;
% vsind_mu = isdf_indices(conj(psir(:,nv-nv_oper+1:nv)), ...
%                               (psir(:,nv-nv_ener+1:nv+nc_ener)), optionsISDF);
% vsgzeta_mu = isdf_kernelg(conj(psir(:,nv-nv_oper+1:nv)), ...
%             (psir(:,nv-nv_ener+1:nv+nc_ener)), vsind_mu, gvec, vol);
[vsind_mu, vsgzeta_mu] = isdf_main('vs', conj(psir(:,nv-nv_oper+1:nv)), ...
        (psir(:,nv-nv_ener+1:nv+nc_ener)), gvec, vol, optionsISDF);
vsgzeta_mu = conj(vsgzeta_mu);
timeforVS = toc(startVS);


startSS = tic;
ssrzeta_mu = zeros(nr, vsrank_mu);
ssgzeta_mu = zeros(ng, ssrank_mu);
optionsISDF.isdfoptions.rank = ssrank_mu;
% [ssrzeta_mu, ssind_mu] = isdf(conj(psir(:,nv-nv_oper+1:nv+nc_oper)), ...
%                               (psir(:,nv-nv_ener+1:nv+nc_ener)), optionsISDF);
% Compute and save aqsch directly here! Prepared for exact CH calculation
% for i = 1:ssrank_mu
%   fftbox1 = reshape(ssrzeta_mu(:, i), gvec.fftgrid);
%   fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * GWinfo.vol;
%   ssgzeta_mu(:, i) = get_from_fftbox(gvec.idxnz, fftbox1, gvec.fftgrid);
% end
% ssgzeta_mu = conj(ssgzeta_mu);
% ssind_mu = isdf_indices(conj(psir(:,nv-nv_oper+1:nv+nc_oper)), ...
%                               (psir(:,nv-nv_ener+1:nv+nc_ener)), optionsISDF);
% ssgzeta_mu = isdf_kernelg(conj(psir(:,nv-nv_oper+1:nv+nc_oper)), ...
%             (psir(:,nv-nv_ener+1:nv+nc_ener)), ssind_mu, gvec, vol);
[ssind_mu, ssgzeta_mu] = isdf_main('ss', conj(psir(:,nv-nv_oper+1:nv+nc_oper)), ...
        (psir(:,nv-nv_ener+1:nv+nc_ener)), gvec, vol, optionsISDF);

ssgzeta_mu = conj(ssgzeta_mu);

timeforSS = toc(startSS);
timeforISDF = toc(startISDF);





sres = zeros(n_ener, 1);
if testflag1
  load I_eps_array.mat
  load coeff_real.mat
  load sres_elements.mat  
end

MWM_occ = zeros(nv_ener, nv_oper, nfreq_real);
MWM_unocc = zeros(nc_ener, nc_oper, nfreq_real);
vcVvc = vcgzeta_mu' * Dcoul * vcgzeta_mu / vol;
vcVnn = vcgzeta_mu' * Dcoul * ssgzeta_mu / vol;

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
      MWM_occ(ibandener_count, ibandoper_Mg, ifreq) ...
      = Mg(:, ibandoper_Mg)' * temp(:, ibandoper_Mg);
    end
  end
  % Calculate for both i and n are unoccupied states.
  % Sae it into MWM_unocc
  for ibandener_count = 1:nc_ener 
    ibandener = bandtocal_unocc(ibandener_count);
    Mg = psir(ssind_mu, ibandener) .* conj(psir(ssind_mu, nv+1:nv+nc_oper));
    Mg = vcVnn * Mg;
    temp = WKernel * Mg;
    for ibandoper_Mg = 1:nc_oper
      MWM_unocc(ibandener_count, ibandoper_Mg, ifreq) ...
      = Mg(:, ibandoper_Mg)' * temp(:, ibandoper_Mg);
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
    x = (ev(ibandener) - ev(ibandoper)) + TOL_SMALL;
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
