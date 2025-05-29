function [ind_mu, zeta_mu] = isdf_main(type, Phig, nlist, mlist, gvec, vol, optionsISDF)
  switch type
    case 'vc'
      fName = 'isdf_typevc.mat';
    case 'vs'
      fName = 'isdf_typevs.mat';
    case 'ss'
      fName = 'isdf_typess.mat';
    otherwise
      error("type should be 'vc', 'vs', or 'ss'");    
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varification of inputs
  flag = isdf_checkInputs(fName, type, Phig, nlist, mlist, gvec, vol, optionsISDF);
% Read or Calculate
  if flag
    % success, read and out
    fprintf("Read and out.\n");
    load(fName, 'ind_mu', 'zeta_mu');
  else
    fprintf("Calculate and out.\n");
    % Calculate wavefunction in real space
    nr = gvec.nfftgridpts; nb = size(Phig, 2); 
    Phig_ = Phig * sqrt(vol);
    Phir = zeros(nr, nb);
    for iband = 1:nb
      fftbox1 = put_into_fftbox(Phig_(:, iband), gvec.idxnz, gvec.fftgrid);
      fftbox1 = nr / vol * do_FFT(fftbox1, gvec.fftgrid, 1);
      Phir(:, iband) = reshape(fftbox1, nr, []);
    end
    % Do isdf
    psi = conj(Phir(:, nlist));
    phi = Phir(:, mlist);
    ind_mu = isdf_indices(psi, phi, optionsISDF);
    zeta_mu = isdf_kernelg(psi, phi, ind_mu, gvec, vol);
    % Failed, do calculation
    ISDFinputs = struct('Phir', Phig, 'nlist', nlist, 'mlist', mlist, 'gvec', gvec, 'vol', vol, 'optionsISDF', optionsISDF);
    save(fName, 'ISDFinputs', 'ind_mu', 'zeta_mu');
  end
end % function
