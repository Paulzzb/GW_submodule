function Ex = gw_x(GWinfor, config, options)
% This file will be final version for the full-frequency GW calculation.
% Generally speaking, user should not directly call this function, instead
% call this indirectly by calling through gwCalculation->gw_fullfreq_cd.

% Description
%   Calculate the exchange part of the self-energies (or frequency-free part).
%   This code only computes the integral on the imaginary axis,  

% Parameters
%   Input:
%   GWinfor: class @GWinfor, contains ground-state information.
%   options: class @GWOptions, contains necessary parameters for the calculation.
%   Output:
%   Ex: the integral on imaginary axis, size nbandener*1.

% Structure 
%   The code is organized as followed:
%   0. Initialize inputs.
%   1. if isdf, do isdf.
%   2. Calculate the self-energies
%

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
    fieldname = nameConstants{i};
    value = default_Constant.(fieldname);    
    if ~isempty(value)
      strEval = sprintf('%s = %.16f;', fieldname, value);
      eval(strEval);
    end
end

ng = GWinfor.gvec.ng; % corresponds to Dcoul

% In the GWinfor, the coulomb matrix is in Rydberg unit.
Dcoul = spdiags(GWinfor.coulG, 0, ng, ng) * ry2ev;
Dcoul(1,1) = GWinfor.coulG0 * ry2ev;

gvec = GWinfor.gvec;
psir = GWinfor.psir;


if (options.isISDF)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Start ISDF
  % Generally, we have relations 
  %   Mr{xx} = xxrzeta_mu * Mrxx_mu;
  startVS = tic;
  optISDF = options;

  
  [vsind_mu, vsgzeta_mu] = isdf_main('vs', psir, nv-nv_oper+1:nv, ...
      nv-nv_ener+1:nv+nc_ener, gvec, vol, optISDF);
  vsgzeta_mu = conj(vsgzeta_mu);
  timeforVS = toc(startVS);
  fprintf('ISDF time = %.4f.\n', timeforVS);


  psirvs = psir(vsind_mu, :);
  Phivs = psirvs(:, nv-nv_oper+1:nv); 
  vsDcoulvs = vsgzeta_mu' * Dcoul * vsgzeta_mu;
  Sigma_x  = vsDcoulvs .* conj(Phivs * Phivs'); 

  Psivs = conj(psirvs(:, nv-nv_ener+1:nv+nc_ener)); 
  Ex  = Psivs' * Sigma_x * Psivs / vol;
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
