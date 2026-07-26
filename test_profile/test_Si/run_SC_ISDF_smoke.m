% Run from test_Si: run_SC_ISDF_smoke
% Smoke test for isdf.SC_ISDF using an adaptive checkpoint on the unit cell.

function run_SC_ISDF_smoke()
  here = pwd;
  case_dir = fileparts(mfilename('fullpath'));
  gw_root = fileparts(fileparts(case_dir));
  service_dir = fullfile(gw_root, 'service');

  cd(service_dir);
  addpath(genpath(service_dir));
  rehash;
  cd(gw_root);
  QPstartup;

  cd(case_dir);
  def = filename_map();
  stage_path = fullfile('SAVE', def.stage);
  if ~isfile(stage_path)
    stage_path = 'test_relay_stage.mat';
  end
  assert(isfile(stage_path), 'Missing relay stage. Run input_driver first.');
  relay.stage_from_db(stage_path);
  relay.restore();

  ck = fullfile('SAVE', 'isdf_adaptive_checkpoint_vc_id1.mat');
  assert(isfile(ck), 'Missing %s', ck);

  uc_fft = int32(FFT.get().fftgrid(:).');
  ratio = int32([1, 1, 1]);
  [id_sc, info] = isdf.SC_ISDF(ck, ratio);
  d = isdf.get(id_sc);
  S = load(ck, 'isdf_data');
  uc = S.isdf_data;

  assert(double(d.nisdf) == double(uc.nisdf), 'nisdf mismatch for 1x1x1 ratio');
  assert(max(abs(d.R_sampling_RLU - uc.R_sampling_RLU), [], 'all') < 1e-9, ...
    'R_sampling mismatch for 1x1x1 ratio');
  assert(numel(d.bundle_struct.fine_grid_lin) == double(d.nisdf), 'fine_grid_lin size');
  fprintf('run_SC_ISDF_smoke: 1x1x1 identity OK (id=%d, nisdf=%d)\n', double(id_sc), double(d.nisdf));

  ratio222 = int32([2, 2, 2]);
  sc_fft = int32(double(uc_fft) .* double(ratio222));
  sc_isdf_smoke_scale_fft(sc_fft);
  out_mat = fullfile('SAVE', 'isdf_sc_vc_2_2_2.mat');
  [id_sc2, info2] = isdf.SC_ISDF(ck, ratio222, 'SavePath', out_mat);
  d2 = isdf.get(id_sc2);
  expected_n = double(uc.nisdf) * 8;
  assert(double(d2.nisdf) == expected_n, 'nisdf supercell mismatch');
  assert(isfile(out_mat), 'SavePath missing');
  fprintf('run_SC_ISDF_smoke: 2x2x2 replication OK (id=%d, nisdf=%d, saved %s)\n', ...
    double(id_sc2), double(d2.nisdf), out_mat);
  fprintf('run_SC_ISDF_smoke: ALL OK\n');
  cd(here);
end

function sc_isdf_smoke_scale_fft(sc_fftgrid)
  fft_data = FFT.get();
  symm_data = symmetry.get();
  nsym = double(symm_data.nsym);
  sc_fftgrid = int32(sc_fftgrid(:)).';
  nr = prod(double(sc_fftgrid));
  fft_m = FFT.base.FFT_m(sc_fftgrid, nsym);
  [I, J, K] = ndgrid(0:sc_fftgrid(1)-1, 0:sc_fftgrid(2)-1, 0:sc_fftgrid(3)-1);
  Rgrid_RLU = double([I(:), J(:), K(:)]);
  fft_m.Rgrid_RLU = Rgrid_RLU;
  nsym_tot = symm_data.nsym;
  is_t_rev = symm_data.is_t_rev;
  for is = 1:int32(nsym_tot + 1)
    if is <= int32(nsym_tot / (1 + is_t_rev))
      mtrx_RLU_R = symm_data.rot_mtrx_RLU_R(:, :, is);
      sign_factor = 1;
    elseif is <= int32(nsym_tot)
      mtrx_RLU_R = symm_data.rot_mtrx_RLU_R(:, :, is);
      sign_factor = -1;
    else
      mtrx_RLU_R = -eye(3);
      sign_factor = 1;
    end
    M2 = double(sign_factor * mtrx_RLU_R);
    M2_r_RLU = (Rgrid_RLU * M2);
    M2_r_RLU = round(M2_r_RLU);
    g = double(sc_fftgrid);
    iv_mod = int32(mod(M2_r_RLU + g, g));
    i4 = 1 + iv_mod(:,1) + iv_mod(:,2)*sc_fftgrid(1) + iv_mod(:,3)*sc_fftgrid(1)*sc_fftgrid(2);
    if is == int32(nsym_tot + 1)
      fft_m.R_rot_inv = int32(i4);
    else
      fft_m.R_rot(:, is) = int32(i4);
    end
  end
  FFT.save2mod(fft_m);
end
