% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function symm_m = driver(data)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % extract from data
  a1a2a3 = data.sys.supercell';
  inv_a1a2a3 = inv(a1a2a3);
  b1b2b3 = 2*pi * inv_a1a2a3';
  fftgrid = [data.sys.n1, data.sys.n2, data.sys.n3];
  %
  syms_in = data.syms;
  is_t_rev = syms_in.is_t_rev;
  nrot = syms_in.nrot;
  nsym = syms_in.nsym;
  mtrx_RLU_G = syms_in.mtrx;
  % init
  symm_m = symmetric.symm_m(nsym, nrot, is_t_rev);
  nsym = symm_m.nsym;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read rotation matrix : RLU on G grid
  for i = 1:int32( single(nsym)/ (1+is_t_rev) )
    symm_m.rot_mtrx_RLU_G(:, :, i) = mtrx_RLU_G{i};
    % single( symm_m.rot_mtrx_RLU_G(:, :, i) );
  end
  if (is_t_rev == 1)
    for i = 1:(nsym/2)
      symm_m.rot_mtrx_RLU_G(:, :, i+nsym/2) = - symm_m.rot_mtrx_RLU_G(:, :, i);
    end
  end

  % Construct rotation in Card .and. RLU on Rgrid
  for irot = 1:nsym
    mtrx_RLU = symm_m.rot_mtrx_RLU_G{irot};
    % Since we assume all vector are row vector, operators are on r.h.s
    mtrx_Cart = a1a2a3 * mtrx_RLU * inv_a1a2a3;
    mtrx_Cart = single ( mtrx_Cart );
    symm_m.rot_mtrx_Cart(:, :, irot) = mtrx_Cart;
    symm_m.rot_mtrx_RLU_R{irot} = (a1a2a3./fftgrid)' * mtrx_Cart * inv(a1a2a3./fftgrid)'; 
  end

  % Construct inverse rotation index
  inv_rot_index = int32( zeros(nsym, 1)-1 );
  for irot = 1:nsym
    inv_rot_Cart_irot = inv(symm_m.rot_mtrx_Cart(:, :, irot));
    for jrot = 1:nsym
      rot_Cart_jrot = symm_m.rot_mtrx_Cart(:, :, jrot);
      if norm(inv_rot_Cart_irot - rot_Cart_jrot) < 1e-5
        inv_rot_index(irot) = jrot;
      end
    end
  end
  symm_m.inv_rot_index = inv_rot_index;

  % put into persistent variable
  symmetry.save2mod(symm_m);
end