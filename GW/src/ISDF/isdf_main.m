function [ind_mu, zeta_mu] = isdf_main(type, psi, phi, gvec, vol, optionsISDF)
  switch type
    case 'vc'
      indfName = 'vcind_mu.mat';
      zetafName = 'vcgzeta_mu.mat';
    case 'vs'
      indfName = 'vsind_mu.mat';
      zetafName = 'vsgzeta_mu.mat';
    case 'ss'
      indfName = 'ssind_mu.mat';
      zetafName = 'ssgzeta_mu.mat';
    otherwise
      warning("type should be 'vc', 'vs', or 'ss'");    
  end

  if (exist(fullfile(pwd, indfName), 'file') == 2)
    ind_mu = load(indfName);
    ind_mu = ind_mu.ind_mu;
  else
    ind_mu = isdf_indices(psi, phi, optionsISDF);
    save(indfName, 'ind_mu');
  end

  if (exist(fullfile(pwd, zetafName), 'file') == 2)
    zeta_mu = load(zetafName);
    zeta_mu = zeta_mu.zeta_mu;
  else
    zeta_mu = isdf_kernelg(psi, phi, ind_mu, gvec, vol);
    save(zetafName, 'zeta_mu')
  end
end % function
