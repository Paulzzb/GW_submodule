function [ind_mu, zeta_mu] = isdf_main(type, Phir, nlist, mlist, gvec, vol, optionsISDF)
  GWlog('\n==================================================\n');
  GWlog('     ISDF Computation     \n');

  def = filename_map();
  switch type
    case 'vc'
      fName = def.isdfvc;
      msg = 'Occupied states with Unoccupied states';
      ratio = optionsISDF.vcrank_ratio;
    case 'vs'
      fName = def.isdfvs;
      msg = 'Occupied states with All states';
      ratio = optionsISDF.vsrank_ratio;
    case 'ss'
      fName = def.isdfss;
      msg = 'All states with All states';
      ratio = optionsISDF.ssrank_ratio;
    otherwise
      error("type should be 'vc', 'vs', or 'ss'");    
  end
  optionsISDF.isdfoptions.rank = ceil(sqrt(length(nlist) * length(mlist)) * ratio);
  msg = sprintf('ISDF Type: %s\n', msg);
  GWlog(msg, 1);
  GWlog('==================================================\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varification of inputs
  GWlog('[Step 1] Verifying inputs...\n', 2);
  flag = isdf_checkInputs(fName, type, Phir, nlist, mlist, gvec, vol, optionsISDF);
% Read or Calculate
  if flag
    % success, read and out
    msg = sprintf('[Step 2] Cached ISDF result found. Loading from %s...\n', fName);
    GWlog(msg, 2);
    load(fName, 'ind_mu', 'zeta_mu');
    msg = sprintf('Loaded interpolative indices and helper functions.\n');
    GWlog(msg, 2);
  else
    msg = sprintf('[Step 2] No valid cache found. Starting ISDF calculation...\n');
    GWlog(msg, 2);
    % Step 3: Prepare wavefunctions
    msg = sprintf('[Step 3] Preparing wavefunctions in real space...\n');
    GWlog(msg, 2);
    psi = conj(Phir(:, nlist));
    phi = Phir(:, mlist);

    % Step 4: Compute interpolation points
    msg = sprintf('[Step 4] Generating interpolation points...\n');
    GWlog(msg, 2);
    ind_mu = isdf_indices(psi, phi, optionsISDF);
    msg = sprintf('Number of interpolation points: %d\n', length(ind_mu));
    GWlog(msg, 1);

    % Step 5: Compute helper functions
    msg = sprintf('[Step 5] Constructing helper functions...\n');
    GWlog(msg, 2);
    zeta_mu = isdf_kernelg(psi, phi, ind_mu, gvec, vol);
    msg = sprintf('Helper functions constructed.\n');
    GWlog(msg, 2);

    % Step 6: Save results
    msg = sprintf('[Step 6] Saving ISDF results to %s...\n', fName);
    GWlog(msg, 2)
    ISDFinputs = struct('Phir', Phir, 'nlist', nlist, 'mlist', mlist, 'gvec', gvec, 'vol', vol, 'optionsISDF', optionsISDF);
    save(fName, 'ISDFinputs', 'ind_mu', 'zeta_mu');
    msg = sprintf('ISDF results saved.\n');
    GWlog(msg, 2);
  end

  msg = sprintf('ISDF processing completed.\n');
  GWlog(msg)
  msg = sprintf('==================================================\n');
  GWlog(msg)
end % function
