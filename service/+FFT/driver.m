% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function driver(data, config)
  % extract from data
  fftgrid = [data.sys.n1, data.sys.n2, data.sys.n3];
  nr = prod(fftgrid);
  %
  % symm_m is constructed in symmetry/driver.m, and is used in FFT/driver.m
  % to construct R_rot and G_rot.
  symm_m_obj = symmetry.get();
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % init
  nsym = double(symm_m_obj.nsym);
  fft_m = FFT.base.FFT_m(fftgrid, nsym);
  
  % Construct Rgrid_RLU: all FFT grid points in RLU (vectorized)
  [I, J, K] = ndgrid(0:fftgrid(1)-1, 0:fftgrid(2)-1, 0:fftgrid(3)-1);
  Rgrid_RLU = double([I(:), J(:), K(:)]);
  fft_m.Rgrid_RLU = Rgrid_RLU;
  
  % Get rotation matrix in RLU space on R grid (FFT grid)
  % These are obtained from symmetry.driver
  
  % Construct R_rot: for each symmetry operation and each FFT grid point,
  % compute where it maps to under rotation
  nsym_tot = symm_m_obj.nsym;
  is_t_rev = symm_m_obj.is_t_rev;
  
  for is = 1:int32(nsym_tot + 1)
    % Select appropriate rotation matrix
    if is <= int32(nsym_tot / (1 + is_t_rev))
      % Standard rotation: positive sign
      % mtrx_RLU_R = symm_m_obj.rot_mtrx_RLU_G(:, :, is);
      mtrx_RLU_R = symm_m_obj.rot_mtrx_RLU_R(:, :, is);
      sign_factor = 1;
    elseif is <= int32(nsym_tot)
      % Time-reversed operation: negative sign
      % mtrx_RLU_R = symm_m_obj.rot_mtrx_RLU_G(:, :, is);
      mtrx_RLU_R = symm_m_obj.rot_mtrx_RLU_R(:, :, is);
      sign_factor = -1;
    else
      % Spatial inversion (is == nsym_tot + 1)
      % Use negative identity for inversion in R grid
      mtrx_RLU_R = -eye(3);
      sign_factor = 1;
    end
    M2 = double( sign_factor * mtrx_RLU_R );
    % for i1=1:3
    %   for i2=1:3
    %     M2(i1,i2) = M2(i1,i2) * fftgrid(i1) / fftgrid(i2);
    %   end
    % end
    % M2 = diag(1./fftgrid) * M2 * diag(fftgrid);

    % Vectorized computation for all grid points
    M2_r_RLU = (Rgrid_RLU * M2);  % nr x 3
    if norm(M2_r_RLU - round(M2_r_RLU)) > 1e-3
      error('Non-integer mapping found in rotation. Check the rotation matrices and FFT grid.');
    end
    M2_r_RLU = (round(M2_r_RLU));
    iv_mod = int32(mod(M2_r_RLU + fftgrid, fftgrid));  % nr x 3
    i4 = 1 + iv_mod(:,1) + iv_mod(:,2)*fftgrid(1) + iv_mod(:,3)*fftgrid(1)*fftgrid(2);  % nr x 1
    
    if is == int32(nsym_tot + 1)
      % Spatial inversion: ir -> i4
      fft_m.R_rot_inv = int32(i4);
    else
      % Standard symmetry: ir -> i4
      fft_m.R_rot(:, is) = int32(i4);
    end
  end 


  % Construct G_table: G(ig)-G(ig0) --> G(G_table(ig, ig0))


  % Save FFT object to manager cache
  FFT.save2mod(fft_m);

  % Report (r-*): owned by this driver
  g = double(fft_m.fftgrid(:).');
  output.msg('nrs', '----------- FFT -----------');
  output.msg('r', ' FFT grid / Nr           :  [%d %d %d] / %d', ...
    int32(g(1)), int32(g(2)), int32(g(3)), int32(fft_m.nr));
end