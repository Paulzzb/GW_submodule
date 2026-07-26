function hVh = cohsex_vcVnn(id_vc, id_nn, nv, nsum, nbmin, nbmax, ikibz, ispin)
%COHSEX_VCVNN  Cross Coulomb matrix <mu_vc | V | nu_nn> in ISDF mu-space (service +isdf).
%
%   hVh = isdf.cohsex_vcVnn(id_vc, id_nn, nv, nsum, nbmin, nbmax, ikibz, ispin)
%
% Algebra matches GW/src_profile/ISDF/isdf_vcVnn.m (FFT + G-space Coulomb),
% but vc / nn interpolation uses coeff_seper from two isdf_m slots instead of
% ISDFDB ind_xga + GWinfo.psir rows.
%
% Band windows (same spirit as legacy vc + ss descriptors):
%   vc block:  valence 1:nv  and conduction nv+1:nsum
%   nn block:  phinn uses 1:nsum, psinn uses nbmin:nbmax (ISDF pair for ss-type C2)

  default_Constant = constant_map();
  nameConstants = fieldnames(default_Constant);
  for ii = 1:numel(nameConstants)
    eval(sprintf('%s = %.16f;', nameConstants{ii}, default_Constant.(nameConstants{ii})));
  end

  vc_data = isdftest.get(id_vc);
  nn_data = isdftest.get(id_nn);
  coeff_vc = vc_data.coeff_seper;
  coeff_nn = nn_data.coeff_seper;

  rkvc = int32(size(coeff_vc, 1));
  rknn = int32(size(coeff_nn, 1));

  fft_data = FFT.manager('get');
  fftgrid = double(fft_data.fftgrid(:).');
  fft_sz = fftgrid;
  idxnz = lattice.manager('r_lat', 'get');
  idxnz = idxnz.idxnz;
  ng = int32(numel(idxnz));

  d_lat_data = lattice.manager('d_lat', 'get');
  vol = double(d_lat_data.DL_vol);

  coul_data = coulomb.get();
  vcoul_q = double(coul_data.vcoul(:, ikibz));
  if ikibz == 1
    vcoul_q(1) = double(coul_data.vcoul0);
  end
  Dcoul = spdiags(vcoul_q(:) * ry2ev, 0, double(ng), double(ng));

  if nv < 1 || nsum <= nv || nbmax < nbmin
    error('isdf:cohsex_vcVnn:Bands', 'Invalid nv=%d nsum=%d nbmin=%d nbmax=%d.', nv, nsum, nbmin, nbmax);
  end
  if nsum > size(coeff_nn, 2) || nbmax > size(coeff_nn, 2)
    error('isdf:cohsex_vcVnn:Nb', 'nsum or nbmax exceeds coeff_nn nb (%d).', size(coeff_nn, 2));
  end

  phivc = double(coeff_vc(:, 1:nv, ikibz, ispin));
  psivc = double(coeff_vc(:, nv + 1:nsum, ikibz, ispin));
  C2 = (phivc * phivc.') .* (psivc * psivc.');
  if condest(C2) < 1e12
    C2invvc = inv(C2);
  else
    C2invvc = pinv(C2);
  end
  clear C2;

  phinn = double(coeff_nn(:, 1:nsum, ikibz, ispin));
  psinn = double(coeff_nn(:, nbmin:nbmax, ikibz, ispin));
  C2 = (phinn * phinn.') .* (psinn * psinn.');
  if condest(C2) < 1e12
    C2invnn = inv(C2);
  else
    C2invnn = pinv(C2);
  end
  clear C2;

  Nblock = 64;
  C1gvc = zeros(double(ng), double(rkvc));
  Psi_v = local_wf_block(ikibz, ispin, 1, nv);
  Psi_c = local_wf_block(ikibz, ispin, nv + 1, nsum);

  for i = 1:Nblock:double(rkvc)
    if i + Nblock < double(rkvc)
      irange = i:i + Nblock - 1;
    else
      irange = i:double(rkvc);
    end
    tmpr = (Psi_v * phivc(irange, :).') .* (Psi_c * psivc(irange, :).');
    for j = 0:length(irange) - 1
      fftbox1 = reshape(tmpr(:, j + 1), fft_sz);
      fftbox1 = do_FFT(fftbox1, fft_sz, 1) * vol;
      C1gvc(:, i + j) = get_from_fftbox(idxnz, fftbox1, fft_sz);
    end
  end

  C1gnn = zeros(double(ng), double(rknn));
  Psi_ns = local_wf_block(ikibz, ispin, 1, nsum);
  Psi_nb = local_wf_block(ikibz, ispin, nbmin, nbmax);

  for i = 1:Nblock:double(rknn)
    if i + Nblock < double(rknn)
      irange = i:i + Nblock - 1;
    else
      irange = i:double(rknn);
    end
    tmpr = (Psi_ns * phinn(irange, :).') .* (Psi_nb * psinn(irange, :).');
    for j = 0:length(irange) - 1
      fftbox1 = reshape(tmpr(:, j + 1), fft_sz);
      fftbox1 = do_FFT(fftbox1, fft_sz, 1) * vol;
      C1gnn(:, i + j) = get_from_fftbox(idxnz, fftbox1, fft_sz);
    end
  end

  C1vcVC1nn = zeros(double(rkvc), double(rknn));
  for i = 1:Nblock:double(rkvc)
    if i + Nblock < double(rkvc)
      irange = i:i + Nblock - 1;
    else
      irange = i:double(rkvc);
    end
    for j = 1:Nblock:double(rknn)
      if j + Nblock < double(rknn)
        jrange = j:j + Nblock - 1;
      else
        jrange = j:double(rknn);
      end
      C1vcVC1nn(irange, jrange) = C1gvc(:, irange)' * Dcoul * C1gnn(:, jrange) / vol;
    end
  end

  hVh = C2invvc * C1vcVC1nn * C2invnn;
end

function M = local_wf_block(ikibz, ispin, ib_lo, ib_hi)
  nr = prod(double(FFT.manager('get').fftgrid));
  nb = ib_hi - ib_lo + 1;
  M = zeros(nr, nb);
  for j = 1:nb
    ib = ib_lo + j - 1;
    isc = int32([ib, ikibz, 1, ispin]);
    M(:, j) = double(wave_functions.WF_apply_symm(isc));
  end
end
