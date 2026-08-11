% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/19 ZZ

function k_ibz = fold_k_to_BZ( type, kpt, b1b2b3 )
  
  ni = 2;
  tol = 1e-5;

  type = lower(type);
  
  if contains(type, 'rlu')
    kpt_RLU = kpt;
    kpt_Cart = kpt_RLU * b1b2b3';
  elseif contains(type, 'cart')
    kpt_Cart = kpt;
    kpt_RLU = kpt_Cart / b1b2b3';
  else
    output.err('type must be either "rlu" or "cart"');
  end

  dist = norm(kpt_Cart, 2);
  k_ibz_RLU = kpt_RLU;
  k_ibz_Cart = kpt_Cart;
  %
  for i1 = -ni:ni
    for i2 = -ni:ni
      for i3 = -ni:ni
        k_i1i2i3_RLU = kpt_RLU - [i1 i2 i3];
        k_i1i2i3_Cart = k_i1i2i3_RLU * b1b2b3';
        if norm( k_i1i2i3_Cart, 2 ) < dist - tol
          k_ibz_Cart = k_i1i2i3_Cart;
          k_ibz_RLU = k_ibz_Cart / b1b2b3';
        end
      end
    end
  end


  if contains(type, 'rlu')
    k_ibz = k_ibz_RLU;
  elseif contains(type, 'cart')
    k_ibz = k_ibz_Cart;
  end



end % function
