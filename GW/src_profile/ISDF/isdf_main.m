function [ind_mu, zeta_mu] = isdf_main(type, Phir, nlist, mlist, gvec, vol, optionsISDF)
  fprintf('\n==================================================\n');
  fprintf('     ISDF Computation     \n');

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
  optionsISDF.isdfoptions.rank = sqrt(length(nlist) * length(mlist)) * ratio;
  fprintf('ISDF Type: %s\n', msg);
  fprintf('==================================================\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varification of inputs
  fprintf('[Step 1] Verifying inputs...\n');
  flag = isdf_checkInputs(fName, type, Phir, nlist, mlist, gvec, vol, optionsISDF);
% Read or Calculate
  if flag
    % success, read and out
    fprintf('[Step 2] Cached ISDF result found. Loading from %s...\n', fName);
    load(fName, 'ind_mu', 'zeta_mu');
    fprintf('Loaded interpolative indices and helper functions.\n');
  else
    fprintf('[Step 2] No valid cache found. Starting ISDF calculation...\n');
    % Step 3: Prepare wavefunctions
    fprintf('[Step 3] Preparing wavefunctions in real space...\n');
    psi = conj(Phir(:, nlist));
    phi = Phir(:, mlist);

    % Step 4: Compute interpolation points
    fprintf('[Step 4] Generating interpolation points...\n');
    ind_mu = isdf_indices(psi, phi, optionsISDF);
    fprintf('Number of interpolation points: %d\n', length(ind_mu));

    % Step 5: Compute helper functions
    fprintf('[Step 5] Constructing helper functions...\n');
    zeta_mu = isdf_kernelg(psi, phi, ind_mu, gvec, vol);
    fprintf('Helper functions constructed.\n');

    % Step 6: Save results
    fprintf('[Step 6] Saving ISDF results to %s...\n', fName);
    ISDFinputs = struct('Phir', Phir, 'nlist', nlist, 'mlist', mlist, 'gvec', gvec, 'vol', vol, 'optionsISDF', optionsISDF);
    save(fName, 'ISDFinputs', 'ind_mu', 'zeta_mu');
    fprintf('ISDF results saved.\n');
  end

  fprintf('ISDF processing completed.\n');
  fprintf('==================================================\n');
end % function
