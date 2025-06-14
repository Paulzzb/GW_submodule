function Ex = gw_x(GWinfor, options)
% This file will be final version for the full-frequency GW calculation.
% Generally speaking, user should not directly call this function, instead
% call this indirectly by calling through gwCalculation->gw_fullfreq_cd.

% Description
%     Calculate the exchange part of the self-energies (or frequency-free part).
%     This code only computes the integral on the imaginary axis,  

% Parameters
%   Input:
%     GWinfor: class @GWinfor, contains ground-state information.
%     options: class @GWOptions, contains necessary parameters for the calculation.
%   Output:
%     Ex: the integral on imaginary axis, size nbandener*1.

% Structure 
%   The code is organized as followed:
%     0. Initialize inputs.
%     1. if isdf, do isdf.
%     2. Calculate the self-energies
%

nameConstants = fieldnames(options.Constant);
for i = 1:numel(nameConstants)
    fieldname = nameConstants{i};
    value = options.Constant.(fieldname);    
    if ~isempty(value)
      strEval = sprintf('%s = %.16f;', fieldname, value);
      eval(strEval);
    end
end

GWinfor.Z = GWinfor.Z * sqrt(GWinfor.vol); % Turn to \int_V \abs{psi(r)}^2 dr = 1.
Dcoul = spdiags(GWinfor.coulG(:,4), 0, ng, ng) * ry2ev;
gvec = GWinfor.gvec;
Dcoul(1,1) = GWinfor.coulG0 * ry2ev;
fprintf('nr = %d, ng = %d, n_oper = %d\n', nr, ng, n_oper);


if (options.ISDFCauchy.isISDF)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Start ISDF
  % Generally, we have relations 
  %     Mr{xx} = xxrzeta_mu * Mrxx_mu;
  startVS = tic;
  optISDF = options.ISDFCauchy;
  vsrank_mu = ceil(sqrt(nv_oper*n_ener)  * optISDF.vsrank_ratio);
  optISDF.isdfoptions.rank = vsrank_mu;

  psir = zeros(nr, nv+nc);
  for iband = 1:nv+nc
    fftbox1 = put_into_fftbox(GWinfor.Z(:, iband), gvec.idxnz, gvec.fftgrid);
    fftbox1 = gvec.nfftgridpts / GWinfor.vol * do_FFT(fftbox1, gvec.fftgrid, 1);
    psir(:, iband) = reshape(fftbox1, gvec.nfftgridpts, []);
  end

  [vsind_mu, vsgzeta_mu] = isdf_main('vs', GWinfor.Z, nv-nv_oper+1:nv, ...
          nv-nv_ener+1:nv+nc_ener, gvec, vol, optISDF);
  vsgzeta_mu = conj(vsgzeta_mu);
  timeforVS = toc(startVS);
  fprintf('ISDF time = %.4f.\n', timeforVS);


  psirvs = psir(vsind_mu, :);
  Phivs = psirvs(:, nv-nv_oper+1:nv); 
  vsDcoulvs = vsgzeta_mu' * Dcoul * vsgzeta_mu;
  Sigma_x    = vsDcoulvs .* conj(Phivs * Phivs'); 

  Psivs = conj(psirvs(:, nv-nv_ener+1:nv+nc_ener)); 
  Ex    = Psivs' * Sigma_x * Psivs / vol;
else
  Ex = zeros(n_ener);
  for ioper = nv-nv_oper+1:nv
    Mgvn = mtxel_sigma(ioper, GWinfor, options.Groundstate, (nv-nv_ener+1:nv+nc_ener));
    Mgvn = conj(Mgvn);
    W1Mgvn = Dcoul * Mgvn;
    Ex = Ex + Mgvn' * W1Mgvn / vol;
  end
end

Ex = - real(diag(Ex));

end % function
