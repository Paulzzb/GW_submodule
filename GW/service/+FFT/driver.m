% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function driver(data)
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
  fft_m = FFT.FFT_m(fftgrid, nsym);
  
  % Construct Rgrid_RLU: all FFT grid points in RLU
  Rgrid_RLU = single( zeros(fftgrid(1)*fftgrid(2)*fftgrid(3), 3) );
  for i = 0:fftgrid(1)-1
    for j = 0:fftgrid(2)-1
      for k = 0:fftgrid(3)-1
        idx = k*fftgrid(1)*fftgrid(2) + j*fftgrid(1) + i + 1;
        Rgrid_RLU(idx, :) = single([i, j, k]);
      end
    end
  end
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
      mtrx_RLU_R = symm_m_obj.rot_mtrx_RLU_R(:, :, is);
      sign_factor = 1;
    elseif is <= int32(nsym_tot)
      % Time-reversed operation: negative sign
      mtrx_RLU_R = symm_m_obj.rot_mtrx_RLU_R(:, :, is);
      sign_factor = -1;
    else
      % Spatial inversion (is == nsym_tot + 1)
      % Use negative identity for inversion in R grid
      mtrx_RLU_R = -eye(3);
      sign_factor = 1;
    end
    M2 = single( sign_factor * mtrx_RLU_R );

    for ir = 1:nr
      M2_r_RLU = round( single( Rgrid_RLU(ir, :) ) * M2 );
      % Apply periodic boundary conditions with modulo (both operands are 1x3 row vectors)
      iv_mod = int32(mod(M2_r_RLU + fftgrid, fftgrid));
      i4 = 1 + iv_mod(1) + iv_mod(2)*fftgrid(1) + iv_mod(3)*fftgrid(1)*fftgrid(2);
      if is == int32(nsym_tot + 1)
        % Spatial inversion: ir -> i4
        fft_m.R_rot_inv(ir) = int32(i4);
      else
        % Standard symmetry: ir -> i4
        fft_m.R_rot(ir, is) = int32(i4);
      end
    end
  end 


  % Construct G_table: G(ig)-G(ig0) --> G(G_table(ig, ig0))


  % Save FFT object to manager cache
  FFT.save2mod(fft_m);
  
end