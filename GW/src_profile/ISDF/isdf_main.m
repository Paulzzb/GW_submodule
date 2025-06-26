function [ind_mu, zeta_mu] = isdf_main(type, Phir, nlist, mlist, gvec, vol, optionsISDF)
  
  cleanup = QPlog_push('ISDF');


  def = filename_map();

  % Determine ISDF type and related settings
  switch type
    case 'vc'
      fName = def.isdfvc;
      outstr = 'Occupied states and Unoccupied states';
      ratio = optionsISDF.vcrank_ratio;
    case 'vs'
      fName = def.isdfvs;
      outstr = 'Occupied states and All states';
      ratio = optionsISDF.vsrank_ratio;
    case 'ss'
      fName = def.isdfss;
      outstr = 'All states and All states';
      ratio = optionsISDF.ssrank_ratio;
    otherwise
      msg = sprintf("ISDF 'type' should be 'vc', 'vs', or 'ss'");
      QPerror(msg);    
  end

  % Estimate rank for ISDF approximation
  optionsISDF.isdfoptions.rank = ceil(sqrt(length(nlist) * length(mlist)) * ratio);

  msg = sprintf('ISDF started');
  QPlog(msg, 0);


  % Step 1: Verify input validity and check for existing result
  QPlog(sprintf(' Verifying inputs...'), 2);
  flag = isdf_checkInputs(fName, type, Phir, nlist, mlist, gvec, vol, optionsISDF);

  if flag
    % Step 2 (cached): Load existing result
    load(fName, 'ind_mu', 'zeta_mu');
    msg = sprintf('Cached ISDF result found. Loading from %s...', fName);
    QPlog(msg, 2);
  else
    % Step 2 (no cache): Start new ISDF computation
    msg = sprintf('No valid cache found. Starting ISDF calculation...');
    QPlog(msg, 2);

    % Step 3: Read wavefunctions in real space

    psi = conj(Phir(:, nlist));
    phi = Phir(:, mlist);

    % Step 3: Compute interpolation points
    msg = sprintf('Generating interpolation points...');
    QPlog(msg, 2);
    ind_mu = isdf_indices(psi, phi, optionsISDF);

    % Step 4: Construct helper functions (zeta_mu)
    msg = sprintf('Constructing helper functions...');
    QPlog(msg, 2);
    zeta_mu = isdf_kernelg(psi, phi, ind_mu, gvec, vol);
    msg = sprintf('Helper functions constructed successfully.');
    QPlog(msg, 2);

    % Step 5: Save result to file
    msg = sprintf('Saving ISDF results to %s...', fName);
    QPlog(msg, 2);
    ISDFinputs = struct('Phir', Phir, 'nlist', nlist, 'mlist', mlist, ...
              'gvec', gvec, 'vol', vol, 'optionsISDF', optionsISDF);
    save(fName, 'ISDFinputs', 'ind_mu', 'zeta_mu');
    msg = sprintf('ISDF results saved successfully.');
    QPlog(msg, 2);
  end

  rank = optionsISDF.isdfoptions.rank;
  msg = sprintf('ISDF completed, Type: %s, Number of helper functions: %d', outstr, rank);
  QPlog(msg, 1);

end % function
