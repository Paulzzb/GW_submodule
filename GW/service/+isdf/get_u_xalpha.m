function u_xalpha = get_u_xalpha(id, isc, iqrot)

  persistent firsttime N_MAX nsym is_t_rev inv_rot_index R_rot sampling2bundle

  if nargin == 1 && (ischar(id) || (isstring(id) && isscalar(id)))
    cmd = lower(string(id));
    if cmd == "reset"
      firsttime = [];
      N_MAX = [];
      inv_rot_index = [];
      R_rot = [];
      sampling2bundle = [];
      u_xalpha = [];
      return;
    end
  end

  if nargin < 2
    error('isdf:get_u_xalpha:Inputs', 'id and isc are required; use get_u_xalpha(''reset'') to clear persistent state.');
  end

  if isempty(N_MAX)
    N_MAX = isdf.isdf_nmax();
    firsttime = true(N_MAX, 1);
    symm_data = symmetry.get();
    % We are facing the same system, so the symmetry properties are the same.
    nsym = symm_data.nsym;
    is_t_rev = symm_data.is_t_rev;
    inv_rot_index = symm_data.inv_rot_index;
    R_rot = cell(N_MAX, 1);
    sampling2bundle = cell(N_MAX, 1);
  end


  if numel(isc) ~= 4
    error('isdf:get_u_xalpha:isc', ...
      'isc must be a length-4 index vector [ib, ikibz, ikrot, ispin]; got %d elements.', numel(isc));
  end
  if any(~isfinite(isc)) || any(isc < 1) || any(isc ~= floor(isc))
    error('isdf:get_u_xalpha:isc', ...
      'isc entries must be finite positive integers; got [%g %g %g %g].', isc(1), isc(2), isc(3), isc(4));
  end

  if nargin < 3 || isempty(iqrot)
    iqrot = 1;
  else
    if ~isscalar(iqrot) || ~isfinite(iqrot) || iqrot < 1 || iqrot ~= floor(iqrot)
      error('isdf:get_u_xalpha:iqrot', ...
        'iqrot must be a positive integer scalar; got: %s.', mat2str(iqrot));
    end
  end

  if firsttime(id)
    isdf_data = isdf.get(id);
    if isempty(isdf_data.bundle_struct) || ~isa(isdf_data.bundle_struct, 'struct')
      error('isdf:get_u_xalpha:bundle', 'ISDF id=%d has no bundle_struct; run rsymm bundle generation first.', int32(id));
    end
    bs = isdf_data.bundle_struct;
    req = {'R_rot_in_bundle', 'sampling2bundle', 'WF_bundle'};
    for k = 1:numel(req)
      if ~isfield(bs, req{k}) || isempty(bs.(req{k}))
        error('isdf:get_u_xalpha:bundle', 'ISDF id=%d bundle_struct missing or empty field ''%s''.', int32(id), req{k});
      end
    end
    R_rot{id} = bs.R_rot_in_bundle;
    sampling2bundle{id} = bs.sampling2bundle;
    firsttime(id) = false;
  end

  isdf_data = isdf.get(id);
  bs = isdf_data.bundle_struct;

  ib    = isc(1);
  ikibz = isc(2);
  ikrot = isc(3);
  ispin = isc(4);

  n_b = size(bs.WF_bundle, 2);
  n_k = size(bs.WF_bundle, 3);
  n_sp = size(bs.WF_bundle, 4);
  if ib > n_b || ikibz > n_k || ispin > n_sp
    error('isdf:get_u_xalpha:bounds', ...
      'WF_bundle is nb=%d, nk=%d, nspin=%d; isc=(%d,%d,_,%d) is out of range.', ...
      n_b, n_k, n_sp, ib, ikibz, ispin);
  end

  % Use bundle_struct
  wf_ikibz_bundle = bs.WF_bundle(:, ib, ikibz, ispin);
  % Apply S_{ikrot} to wf_ikibz_bundle
  % Since now S_{ikrot} contains sign information, we need to apply it carefully.
  if ikrot ~= 1
    isconj = false;
    if ikrot > nsym / (1 + is_t_rev)
      ikrot = ikrot - nsym / (1 + is_t_rev);
      isconj = true;
    end
    invikrot = inv_rot_index(ikrot);

    wf_ikbz_bundle = wf_ikibz_bundle(R_rot{id}(:, invikrot));
    if isconj
      wf_ikbz_bundle = conj(wf_ikbz_bundle);
    end
  else
    wf_ikbz_bundle = wf_ikibz_bundle;
  end

  % Then apply S_{iqrot} to wf_ikbz_bundle
  if iqrot ~= 1
    wf_ikbz_bundle = wf_ikbz_bundle(R_rot{id}(:, iqrot));
  end
  u_xalpha = wf_ikbz_bundle(sampling2bundle{id});

  % Optional consistency checks under isdf.debug switches.
  if isdf.debug.on('u_xalpha_sampling')
    fft_data = FFT.manager('get');
    symm_data = symmetry.manager('get');
    N_coarse = isdf_data.bundle_struct.N_coarse;
    N_sampling = isdf_data.bundle_struct.N_sampling;
    fftgrid = single(fft_data.fftgrid);
    R_sampling = isdf_data.bundle_struct.R_grid_bundle(isdf_data.bundle_struct.sampling2bundle, :);
    R_sampling2 = isdf_data.R_sampling_RLU;
    isdf.debug.check('u_xalpha_sampling', ...
      @() norm(R_sampling - R_sampling2) > 1e-3, ...
      'The sampling rotation on the bundle is not correct', ...
      'u_xalpha_sampling_rotation');

    Sq_R_sampling = R_sampling * symm_data.rot_mtrx_RLU_R(:, :, iqrot);
    Sq_R_sampling = mod(round(Sq_R_sampling) + fftgrid, fftgrid);
    ind_s2f = 1 + Sq_R_sampling(:, 1) + fftgrid(1) * Sq_R_sampling(:, 2) ...
                       + fftgrid(1) * fftgrid(2) * Sq_R_sampling(:, 3);
    ind_s2f = int32(round(ind_s2f));
    wf = wave_functions.WF_apply_symm(isc);
    wf_Sq_R_sampling = wf(ind_s2f);
    isdf.debug.check('u_xalpha_sampling', ...
      @() norm(u_xalpha(N_coarse + 1:N_sampling) - wf_Sq_R_sampling(N_coarse + 1:N_sampling)) > 1e-3, ...
      'The wavefunction on the sampling is not correct', ...
      'u_xalpha_sampling_wavefunction');
  end

  % Do a test here
  % u_xalpha_test = zeros(size(u_xalpha));
  % %
  % N_coarse = isdf_data.bundle_struct.N_coarse;
  % N_sampling = isdf_data.bundle_struct.N_sampling;
  % fft_data = FFT.manager('get');
  % fftgrid = single(fft_data.fftgrid);
  % symm_data = symmetry.get();
  % %
  % u_ikibz_c = isdf_data.coeff_seper(:, ib, ikibz, ispin);
  % u_ikbz_c = isdf.coeff.isdf_apply_symm_on_coarse(id, u_ikibz_c, ikrot);
  % ind_Sq_R_coarse = isdf_data.R_rot_coarse(:, iqrot);
  % u_xalpha_test(1:N_coarse) = u_ikbz_c(ind_Sq_R_coarse);
  % %
  % u_ikbz = wave_functions.WF_apply_symm(isc);
  % R_f_sampling = isdf_data.bundle_struct.R_grid_bundle(isdf_data.bundle_struct.sampling2bundle(N_coarse+1:N_sampling), :);
  % Sq_R_f_sampling = R_f_sampling * symm_data.rot_mtrx_RLU_R(:, :, iqrot);
  % if norm(Sq_R_f_sampling - round(Sq_R_f_sampling)) > 1e-3
  %   error('The sampling rotation on the bundle is not correct');
  % end
  % Sq_R_f_sampling = single( round(Sq_R_f_sampling) );
  % Sq_R_f_sampling = mod(round(Sq_R_f_sampling) + fftgrid, fftgrid);
  % ind = 1 + Sq_R_f_sampling(:, 1) + fftgrid(1) * Sq_R_f_sampling(:, 2) ...
  %        + fftgrid(1) * fftgrid(2) * Sq_R_f_sampling(:, 3);
  % u_xalpha_test(N_coarse+1:N_sampling) = u_ikbz(ind);
  % if isconj
  %   u_xalpha_test = conj(u_xalpha_test);
  % end
  % if norm(u_xalpha - u_xalpha_test) > 8e-5
  %   error('The wavefunction on the coarse grid is not correct');
  % end
end
