% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function RL_shell_construct(ecut)
% 1. Select G-vectors within the cutoff Ecut, such that
%    |G|^2 < Ecut.
%    This will be saved in r_lattice_m.G_vec_RLU and r_lattice_m.G_vec_Cart.
% 2. Construct shell structure: group G-vectors into shells based on their |G|^2 values.
%    This will be saved in n_g_shell, num_index_in_each_shell,
%    and first_index_of_each_shell in `r_lattice_m'. 
  symm_m = symmetry.manager('get');
  r_lat_m = lattice.manager('r_lat', 'get');
  fft_m = FFT.manager('get');
  %
  fftgrid = fft_m.fftgrid;
  fftgrid_s = single(fftgrid);
  b1b2b3 = r_lat_m.b1b2b3;
  %
  n1 = fftgrid(1);
  n2 = fftgrid(2);
  n3 = fftgrid(3); 

  n1_s = single(n1);
  n2_s = single(n2);
  n3_s = single(n3);
  [gkxind, gkyind, gkzind] = ndgrid( ...
  (0:n1_s-1) - ( (0:n1_s-1) >= n1_s/2 )*n1_s, ...
  (0:n2_s-1) - ( (0:n2_s-1) >= n2_s/2 )*n2_s, ...
  (0:n3_s-1) - ( (0:n3_s-1) >= n3_s/2 )*n3_s);
  gkxind = int32( gkxind(:) );
  gkyind = int32( gkyind(:) );
  gkzind = int32( gkzind(:) );
  
  % Calculate |G|^2 for all G-vectors in the full FFT grid,
  % and select those within the cutoff 
  gkind = [gkxind, gkyind, gkzind];
  gkvec = single(gkind) * b1b2b3';
  gkabs2 = sum(gkvec.^2, 2);
  % Sort according to |G|.^2
  [gkabs2, sort_ind] = sort(gkabs2);
  gkind = gkind(sort_ind, :);
  gkxind = gkxind(sort_ind);
  gkyind = gkyind(sort_ind);
  gkzind = gkzind(sort_ind);

  % ecut =   59.7353365;
  idxnz = int32( find(gkabs2 <= ecut) );
  r_lat_m.ng = length(idxnz);
  gkxind = gkxind(idxnz);
  gkyind = gkyind(idxnz);
  gkzind = gkzind(idxnz);
  gkabs2 = gkabs2(idxnz);

  % gkabs2 = gkabs2(sort_ind);
  r_lat_m.Ggrid_RLU = [gkxind, gkyind, gkzind];
  r_lat_m.Ggrid_Cart = single(r_lat_m.Ggrid_RLU) * b1b2b3';
  % r_lat_m.idxnz = idxnz(sort_ind);

  % Construct shell structure
  first_index_in_each_shell = zeros(r_lat_m.ng, 1);
  num_shell = 0;
  num_index_in_each_shell = zeros(r_lat_m.ng, 1);
  % ig=1
  num_shell = num_shell + 1;
  num_index_in_each_shell(num_shell) = 1;
  first_index_in_each_shell(num_shell) = 1;
  %
  for ig = 2:r_lat_m.ng
    old = gkabs2(ig-1);
    new = gkabs2(ig);
    if abs(new-old) > r_lat_m.tol
      num_shell = num_shell + 1;
      first_index_in_each_shell(num_shell) = ig;
    end
    num_index_in_each_shell(num_shell) = num_index_in_each_shell(num_shell) + 1;
  end
  r_lat_m.n_g_shell = num_shell;
  r_lat_m.num_index_in_each_shell = num_index_in_each_shell(1:num_shell);
  r_lat_m.first_index_in_each_shell = first_index_in_each_shell(1:num_shell);

  % Construct G_rot
  r_lat_m.G_rot = int32( zeros(r_lat_m.ng, symm_m.nsym) );
  for irot = 1:symm_m.nsym 
    mtrx_RLU = symm_m.rot_mtrx_RLU_G(:, :, irot);
    rot_Ggrid_RLU = int32( single( r_lat_m.Ggrid_RLU(:, :) ) * mtrx_RLU );
    for ishell = 1:r_lat_m.n_g_shell
      ig_start = r_lat_m.first_index_in_each_shell(ishell);
      ig_len = r_lat_m.num_index_in_each_shell(ishell);
      ig_end = ig_start + ig_len - 1;
      ig_range = ig_start:ig_end;
      for ig1 = ig_range
        for ig2 = ig_range
          gdiff = rot_Ggrid_RLU(ig1, :) - r_lat_m.Ggrid_RLU(ig2, :); 
          if all( mod(gdiff + fftgrid, fftgrid) == 0 ) 
            r_lat_m.G_rot(ig1, irot) = ig2;
            break
          end
          if ig2 == ig_end
            error('Cannot find the rotated G-vector in the same shell. This should not happen if the tolerance is set properly.');
          end
        end
      end
    end
  end
  
  lattice.manager( 'r_lat', 'save2mod', r_lat_m );

end % EOF