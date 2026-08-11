% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function driver(data, config)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % extract from data
  a1a2a3 = double( data.sys.supercell' );
  inv_a1a2a3 = inv(a1a2a3);
  b1b2b3 = 2*pi * inv_a1a2a3';
  fftgrid = int32( [data.sys.n1, data.sys.n2, data.sys.n3] );
  %
  syms_in = data.syms;
  is_t_rev = syms_in.is_t_rev;
  nrot = syms_in.nrot;
  nsym = syms_in.nsym;
  mtrx_RLU_G = syms_in.mtrx;
  % init
  symm_data = symmetry.base.symm_m(nsym, nrot, is_t_rev);
  nsym = symm_data.nsym;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read rotation matrix : RLU on G grid
  for i = 1: nsym / (1+is_t_rev)
    symm_data.rot_mtrx_RLU_G(:, :, i) = mtrx_RLU_G{i};
    % double( symm_data.rot_mtrx_RLU_G(:, :, i) );
  end
  if (is_t_rev == 1)
    for i = 1:( nsym/(1+is_t_rev) )
      symm_data.rot_mtrx_RLU_G(:, :, i+nsym/2) = - symm_data.rot_mtrx_RLU_G(:, :, i);
    end
  end

  % Construct rotation in Card .and. RLU on Rgrid
  % a1a2a3_scal = a1a2a3 ./ double(fftgrid);
  R_scal = diag(1./double(fftgrid));
  % inv_a1a2a3_scal = inv(a1a2a3_scal);
  for irot = 1:nsym
    mtrx_RLU = double( symm_data.rot_mtrx_RLU_G(:, :, irot) );
    % Since we assume all vector are row vector, operators are on r.h.s
    mtrx_Cart = a1a2a3 * mtrx_RLU * inv_a1a2a3;
    mtrx_Cart = double( mtrx_Cart );
    symm_data.rot_mtrx_Cart(:, :, irot) = mtrx_Cart;
    symm_data.rot_mtrx_RLU_R(:, :, irot) = ...
      R_scal * (a1a2a3)' * mtrx_Cart * inv_a1a2a3' * diag(double(fftgrid));
    % symm_data.rot_mtrx_RLU_R(:, :, irot) =  mtrx_RLU; 
  end

  % Construct inverse rotation index
  inv_rot_index = int32( zeros(nsym, 1)-1 );
  for irot = 1:nsym
    inv_rot_Cart_irot = inv(symm_data.rot_mtrx_Cart(:, :, irot));
    for jrot = 1:nsym
      rot_Cart_jrot = symm_data.rot_mtrx_Cart(:, :, jrot);
      if norm(inv_rot_Cart_irot - rot_Cart_jrot) < 1e-5
        inv_rot_index(irot) = jrot;
      end
    end
  end
  symm_data.inv_rot_index = inv_rot_index;

  % put into persistent variable
  symmetry.save2mod(symm_data);

  % Report (r-*): owned by this driver
  output.msg('nrs', '----------- Symmetry -----------');
  output.msg('r', ' nsym / nrot / t_rev     :  %d / %d / %d', ...
    int32(symm_data.nsym), int32(symm_data.nrot), int32(symm_data.is_t_rev));
end