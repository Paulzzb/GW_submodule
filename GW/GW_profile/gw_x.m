function Ex = gw_x(GWinfor, config)
% gw_x -> Compute the static exchange (Hartree-Fock) self-energy correction.
%
% This function computes the exchange contribution (Σ_x) to the GW
% self-energy, based on either standard formulation or the Interpolative 
% Separable Density Fitting (ISDF) approximation.
%
% Inputs:
%   GWinfor - Struct of class @GWinfor, containing ground-state data:
%       .gvec      : Reciprocal grid info (FFT grid, indices)
%       .vol       : Unit cell volume
%       .coulG     : G-space Coulomb potential
%       .coulG0    : Special treatment of G=0 component
%       .psir      : Wavefunctions in real space
%       .occupation: Orbital occupation numbers
%   config - Struct of class @GWOptions, containing:
%       .SYSTEM.energy_band_index_min / max: band index range
%       .ISDF.isisdf: flag to use ISDF
%       .ISDFCauchy : ISDF options (e.g., rank, seed, method)
%
% Output:
%   Ex - Column vector (n × 1) of exchange corrections for bands n in [min, max]
%
% Note:
%   The exchange energy is computed as:
%     Ex_n = -⟨ψ_n | Σ_x | ψ_n⟩
%   Units: electron-Volts (eV)


msg = sprintf('[Exchange] Start computing Σ_x (exchange part)...\n');
QPlog(msg, 0);
tStart = tic;

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end


% Basic setup
ng = GWinfor.gvec.ng; % corresponds to Dcoul
nbmin = config.SYSTEM.energy_band_index_min;
nbmax = config.SYSTEM.energy_band_index_max;
nv = find(GWinfor.occupation > 1 - TOL_SMALL, 1, 'last');
vol = GWinfor.vol;
gvec = GWinfor.gvec;
psir = GWinfor.psir;

msg = sprintf('[Exchange] Using band range [%d, %d], %d valence bands detected.\n', ...
         nbmin, nbmax, nv);
QPlog(msg);

% Construct Coulomb matrix in sparse form, unit: eV
Dcoul = spdiags(GWinfor.coulG, 0, ng, ng) * ry2ev;
Dcoul(1,1) = GWinfor.coulG0 * ry2ev;



if (config.ISDF.isisdf)
  msg = sprintf('[Exchange] Using ISDF approximation for Σ_x.\n');
  QPlog(msg, 0);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tISDF = tic;
  
  optISDF = config.ISDFCauchy;
  
  % Perform ISDF
  [vsind_mu, vsgzeta_mu] = isdf_main('vs', psir, 1:nv, ...
      nbmin:nbmax, gvec, vol, optISDF);
  vsgzeta_mu = conj(vsgzeta_mu);
  msg = sprintf('[Exchange] ISDF done in %.2f seconds.\n', toc(tISDF));
  QPlog(msg, 1);


  % Compute Σ_x
  psirvs = psir(vsind_mu, :);
  Phivs = psirvs(:, 1:nv); 
  vsDcoulvs = vsgzeta_mu' * Dcoul * vsgzeta_mu;
  Sigma_x  = vsDcoulvs .* conj(Phivs * Phivs'); 

  Psivs = conj(psirvs(:, nbmin:nbmax)); 
  Ex  = Psivs' * Sigma_x * Psivs / vol;

else
% --- Standard exchange calculation without ISDF ---
  msg = sprintf('[Exchange] Using standard Σ_x calculation.\n');
  QPlog(msg, 0);
  Ex = zeros(nbmax-nbmin+1);
  tStandard = tic;
  for ioper = 1:nv
    Mgvn = mtxel_sigma(ioper, GWinfor, nbmin:nbmax);
    Mgvn = conj(Mgvn);
    W1Mgvn = Dcoul * Mgvn;
    Ex = Ex + Mgvn' * W1Mgvn / vol;
  end
  msg = sprintf('[Exchange] Standard loop completed in %.2f seconds.\n', toc(tStandard));
  QPlog(msg, 1);
end

Ex = - real(diag(Ex));
msg = sprintf('[Exchange] Finished. Total time: %.2f seconds.\n', toc(tStart));
QPlog(msg, 0);

end % function
