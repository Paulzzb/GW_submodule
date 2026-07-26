function k_ibz = bz2ibz( type, kpt, b1b2b3 )
  
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
    error('type must be either "rlu" or "cart"');
  end
  % k_ibz = kpt_Cart;
  % k_ibz(:,1) = mod(kpt_Cart(:,1),1);
  % k_ibz(:,2) = mod(kpt_Cart(:,2),1);
  % k_ibz(:,3) = mod(kpt_Cart(:,3),1);
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
        % k_ibz = [k_ibz; kpt_Cart + i1*b1b2b3(1,:) + i2*b1b2b3(2,:) + i3*b1b2b3(3,:)];
        end
      end
    end
  end

  % mask = abs(abs(k_ibz_RLU) - 0.5) < tol | abs(k_ibz_RLU) < tol;
  % if ~all(mask)
  %   return
  % end
  % idx = find(abs(abs(k_ibz_RLU) - 0.5) < tol, 1);
  % % idx = find(abs(k_ibz_Cart) > tol, 1);
  % if ~isempty(idx)
  %   if k_ibz_RLU(idx) > tol 
  %   % if k_ibz_Cart(idx) > tol 
  %     k_ibz_RLU = -k_ibz_RLU;
  %     k_ibz_Cart = -k_ibz_Cart;
  %   end
  % end

  if contains(type, 'rlu')
    k_ibz = k_ibz_RLU;
  elseif contains(type, 'cart')
    k_ibz = k_ibz_Cart;
  end



end % function
