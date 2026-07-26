function gw = demo_GW_scaling(out)
%DEMO_GW_SCALING  Formal Gamma COHSEX from demo_isdf_scaling report (no main-path GW).
%
%   gw = demo_GW_scaling(out)
%
%   Algebra follows GW/GW_profile/gw_cohsex.m (ISDF block, ~L91–183).
%   Uses out.isdf.vcVvc, out.isdf.vcVnn, psixga (wf.c at mu points) from demo_isdf_scaling.

  req_out = {'data_sc', 'config_sc', 'isdf', 'UC_save_dir', 'k1k2k3', 'data_uc'};
  for k = 1:numel(req_out)
    if ~isfield(out, req_out{k})
      error('demo_GW_scaling:field', 'Missing field out.%s.', req_out{k});
    end
  end

  default_Constant = constant_map();
  for nm = fieldnames(default_Constant).'
    eval(sprintf('%s = %.16f;', nm{1}, default_Constant.(nm{1})));
  end

  data_sc = out.data_sc;
  config = out.config_sc;
  isdf_rep = out.isdf;
  if ~isfield(isdf_rep, 'vc') || ~isfield(isdf_rep, 'nn')
    error('demo_GW_scaling:isdf', 'out.isdf must contain .vc and .nn.');
  end
  rep_vc = isdf_rep.vc;
  rep_nn = isdf_rep.nn;
  req_case = {'helperqG', 'CCHq', 'vcoul_q', 'psixga'};
  for k = 1:numel(req_case)
    if ~isfield(rep_vc, req_case{k}) || ~isfield(rep_nn, req_case{k})
      error('demo_GW_scaling:case', 'out.isdf.vc/nn missing %s; rerun demo_isdf_scaling.', req_case{k});
    end
  end

  nbmin = config.SYSTEM.energy_band_index_min;
  nbmax = config.SYSTEM.energy_band_index_max;
  nb = size(data_sc.psig{1}, 2);
  nsum = min(config.SYSTEM.number_bands_in_summation, nb);

  occ = data_sc.occupation;
  if ndims(occ) == 2
    occ = reshape(occ, size(occ, 1), 1, 1);
  end
  nv = find(occ(:, 1, 1) > 1 - TOL_SMALL, 1, 'last');
  if isempty(nv)
    error('demo_GW_scaling:nv', 'Cannot determine nv from data_sc.occupation.');
  end
  if nsum <= nv
    error('demo_GW_scaling:nsum', 'Need nsum > nv (got nv=%d, nsum=%d).', nv, nsum);
  end

  ev = squeeze(double(data_sc.ev(:, 1, 1)));

  if isfield(isdf_rep, 'vcVvc') && ~isempty(isdf_rep.vcVvc)
    vcVvc = double(isdf_rep.vcVvc);
  elseif isfield(rep_vc, 'tildeVq') && ~isempty(rep_vc.tildeVq)
    vcVvc = double(rep_vc.tildeVq);
  else
    error('demo_GW_scaling:vcVvc', 'Missing out.isdf.vcVvc; rerun demo_isdf_scaling.');
  end
  if ndims(vcVvc) == 3
    vcVvc = vcVvc(:, :, 1);
  end

  if isfield(isdf_rep, 'vcVnn') && ~isempty(isdf_rep.vcVnn)
    vcVnn = double(isdf_rep.vcVnn);
  else
    error('demo_GW_scaling:vcVnn', 'Missing out.isdf.vcVnn; rerun demo_isdf_scaling.');
  end

  Nmu_vc = size(rep_vc.helperqG, 2);
  Nmu_nn = size(rep_nn.helperqG, 2);
  if size(rep_vc.psixga, 1) ~= Nmu_vc || size(rep_nn.psixga, 1) ~= Nmu_nn
    error('demo_GW_scaling:psixga', 'psixga row count != Nmu; rerun demo_isdf_scaling.');
  end
  if nsum > size(rep_vc.psixga, 2) || nsum > size(rep_nn.psixga, 2)
    error('demo_GW_scaling:psixga', 'psixga missing bands for nsum=%d.', nsum);
  end
  if nbmax > size(rep_nn.psixga, 2)
    error('demo_GW_scaling:psixga', 'psixga missing bands for nbmax=%d.', nbmax);
  end

  fprintf(['\n[demo_GW_scaling] formal COHSEX | nv=%d nsum=%d bands [%d,%d] ', ...
    'Nmu_vc=%d Nmu_nn=%d | COmegaC=%s\n'], nv, nsum, nbmin, nbmax, Nmu_vc, Nmu_nn, ...
    local_cohsex_comega_label(config));

  psixgavc = rep_vc.psixga(:, 1:nsum);
  psixgass = rep_nn.psixga;

  % COmegaCresult  (gw_cohsex.m L117–140)
  startinveps = tic;
  if local_use_cauchy(config)
    Phi = psixgavc(:, 1:nv);
    Psi = psixgavc(:, nv + 1:nsum);
    evOcc = ev(1:nv);
    evUnocc = ev(nv + 1:nsum);
    [COmegaCresult, ~, ~] = COmegaCstar(Phi, Psi, evOcc, evUnocc, ...
      local_cauchy_options(config));
    COmegaCresult = COmegaCresult * 4;
  else
    COmegaCresult = zeros(Nmu_vc, Nmu_vc);
    scal = 4.0;
    for ind_nv = 1:nv
      Mgvc = conj(psixgavc(:, ind_nv)) .* psixgavc(:, nv + 1:nsum);
      Mgvc = conj(Mgvc);
      eden = 1 ./ (ev(ind_nv) - ev(nv + 1:nsum));
      COmegaCresult = COmegaCresult + scal * Mgvc * diag(eden) * Mgvc';
    end
  end
  epsg_main = inv(COmegaCresult) - vcVvc;
  timeforinveps = toc(startinveps);
  fprintf('[demo_GW_scaling] Time for inveps = %.4f s\n', timeforinveps);

  % Sigma operators (gw_cohsex.m L148–163)
  startSigma = tic;
  Phivs = psixgass(:, 1:nv);
  Phiss = psixgass(:, 1:nsum);
  epsvcDcoulss = epsg_main \ vcVnn;
  W1_mu_1 = vcVnn' * epsvcDcoulss;
  Sigma_sex_x = W1_mu_1 .* (Phivs * Phivs');
  Sigma_coh = W1_mu_1 .* (Phiss * Phiss');
  timeForSigma = toc(startSigma);
  fprintf('[demo_GW_scaling] Sigma operator time = %.4f s\n', timeForSigma);

  % Self-energies (gw_cohsex.m L173–178)
  startSelfE = tic;
  Psivs = conj(psixgass(:, nbmin:nbmax));
  Psiss = conj(psixgass(:, nbmin:nbmax));
  Esx_x = -Psivs' * Sigma_sex_x * Psivs;
  Ecoh = 0.5 * Psiss' * Sigma_coh * Psiss;
  timeforSelfE = toc(startSelfE);
  fprintf('[demo_GW_scaling] Self-energies time = %.4f s\n', timeforSelfE);

  gw = struct();
  gw.Esx_x = real(diag(Esx_x));
  gw.Ecoh = real(diag(Ecoh));
  gw.wall_s = timeforinveps + timeForSigma + timeforSelfE;
  gw.nv = nv;
  gw.nsum = nsum;
  gw.Nmu_vc = Nmu_vc;
  gw.Nmu_nn = Nmu_nn;
end

function tf = local_use_cauchy(config)
  tf = false;
  if isfield(config, 'ISDFCauchy') && isstruct(config.ISDFCauchy) ...
      && isfield(config.ISDFCauchy, 'isCauchy')
    tf = logical(config.ISDFCauchy.isCauchy);
  end
end

function label = local_cohsex_comega_label(config)
  if local_use_cauchy(config)
    label = 'cauchy';
  else
    label = 'direct';
  end
end

function opt = local_cauchy_options(config)
  opt = struct('froErr', 1e-6, 'MaxIter', 10);
  if ~isfield(config, 'ISDFCauchy') || ~isstruct(config.ISDFCauchy)
    return;
  end
  ic = config.ISDFCauchy;
  if isfield(ic, 'froErr') && ~isempty(ic.froErr)
    opt.froErr = ic.froErr;
  end
  if isfield(ic, 'MaxIter') && ~isempty(ic.MaxIter)
    opt.MaxIter = ic.MaxIter;
  end
  if isfield(ic, 'optionsCauchy') && isstruct(ic.optionsCauchy) && ~isempty(ic.optionsCauchy)
    oc = ic.optionsCauchy;
    if isfield(oc, 'froErr') && ~isempty(oc.froErr)
      opt.froErr = oc.froErr;
    end
    if isfield(oc, 'MaxIter') && ~isempty(oc.MaxIter)
      opt.MaxIter = oc.MaxIter;
    end
  end
end
