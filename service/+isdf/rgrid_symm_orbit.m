function [N_r_orbit, Nr, Rgrid_RLU, Rgrid_rep_RLU, ir2rep, ir2rot] = rgrid_symm_orbit()
%RGRID_SYMM_ORBIT Build symmetry-orbit structure on FFT real-space grid.
%
% Output convention:
%   N_r_orbit : number of orbits on R-grid
%   Nr        : total number of R-grid points
%   Rgrid_RLU : full R-grid points in RLU, size [Nr, 3]
%   Rgrid_rep_RLU : orbit representatives in RLU, size [N_r_orbit, 3]
%   ir2rep    : representative index for each ir, size [Nr, 1]
%   ir2rot    : rotation index for each ir, size [Nr, 1]
%
% Reconstruction relation:
%   Rgrid_RLU(ir, :) = Rgrid_rep_RLU(ir2rep(ir), :) * rot_mtrx_RLU_R(:, :, ir2rot(ir))

  fft_data = FFT.get();
  symm_data = symmetry.get();

  Rgrid_RLU = int32(fft_data.Rgrid_RLU);
  Nr = int32(size(Rgrid_RLU, 1));
  fftgrid = int32(fft_data.fftgrid);
  rot_mtrx_RLU_R = double(symm_data.rot_mtrx_RLU_R);
  nsym = int32(size(rot_mtrx_RLU_R, 3));
  is_t_rev = symm_data.is_t_rev;

  if nsym <= 0
    error('rgrid_symm_orbit:NoSymmetry', 'symmetry.rot_mtrx_RLU_R is empty.');
  end

  ir2rep = int32(zeros(Nr, 1));
  ir2rot = int32(zeros(Nr, 1));
  Rgrid_rep_RLU = int32(zeros(0, 3));

  rep_count = int32(0);
  tol = 1e-3;

  for ir = 1:Nr
    if ir2rep(ir) ~= 0
      continue;
    end

    rep_count = rep_count + 1;
    r_rep = double(Rgrid_RLU(ir, :));
    Rgrid_rep_RLU(rep_count, :) = int32(r_rep); %#ok<AGROW>

    for irot = 1 : nsym / (1 + is_t_rev)
      r_img = r_rep * rot_mtrx_RLU_R(:, :, irot);
      r_img_round = round(r_img);
      if norm(r_img - r_img_round) > tol
        error('rgrid_symm_orbit:NonIntegerRotate', ...
          'R-grid rotation gives non-integer image (rep=%d, irot=%d).', ...
          rep_count, irot);
      end

      r_mod = mod(int32(r_img_round) + fftgrid, fftgrid);
      ir_img = int32(1 + r_mod(1) + r_mod(2) * fftgrid(1) + ...
        r_mod(3) * fftgrid(1) * fftgrid(2));

      if ir2rep(ir_img) == 0
        ir2rep(ir_img) = rep_count;
        ir2rot(ir_img) = irot;
      end
    end
  end

  if any(ir2rep == 0) || any(ir2rot == 0)
    error('rgrid_symm_orbit:IncompleteCover', ...
      'Some R-grid points are not covered by symmetry orbits.');
  end

  N_r_orbit = rep_count;

  % Consistency check against the required reconstruction relation.
  for ir = 1:Nr
    irep = ir2rep(ir);
    irot = ir2rot(ir);
    r_rhs = double(Rgrid_rep_RLU(irep, :)) * rot_mtrx_RLU_R(:, :, irot);
    r_rhs = mod(int32(round(r_rhs)) + fftgrid, fftgrid);
    if any(r_rhs ~= Rgrid_RLU(ir, :))
      error('rgrid_symm_orbit:RelationMismatch', ...
        'Orbit relation mismatch at ir=%d.', ir);
    end
  end
end
