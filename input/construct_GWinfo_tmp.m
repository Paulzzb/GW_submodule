function GWinfo = construct_GWinfo_tmp(GWinfo)
  
  tmp_devel = struct();

  gvec = GWinfo.gvec;
  bz_samp = GWinfo.bz_samp;
  
  % Real lattice vectors and reciprocal lattice vectors

  a1a2a3 = double( GWinfo.supercell' );
  b1b2b3 = double( 2*pi * inv(a1a2a3)');

  % Ggrid, Rgrid
  ng = gvec.ng;
  fftgrid = gvec.fftgrid;
  nr = gvec.nfftgridpts;
  %
  Ggrid_RLU = round(gvec.components);
  Ggrid_Cart = double(Ggrid_RLU) * b1b2b3';
  %
  Rgrid_RLU = double( zeros(nr, 3) );
  count = 0;
  for i3 = 0:fftgrid(3)-1
    for i2 = 0:fftgrid(2)-1
      for i1 = 0:fftgrid(1)-1
        count = count + 1;
        Rgrid_RLU(count, :) = [i1, i2, i3];
      end
    end
  end
  Rgrid_Cart = Rgrid_RLU * a1a2a3';
  tmp_devel.Rgrid_Cart = Rgrid_Cart;
  tmp_devel.Ggrid_Cart = Ggrid_Cart;
  tmp_devel.Ggrid_RLU = Ggrid_RLU;
  tmp_devel.Rgrid_RLU = Rgrid_RLU;
  tmp_devel.a1a2a3 = a1a2a3;
  tmp_devel.b1b2b3 = b1b2b3;
  tmp_devel.fftgrid = fftgrid;
  tmp_devel.ng = ng;
  tmp_devel.nr = nr;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % rotation group
  % Rotation matrices
  nsym = GWinfo.symminfo.nsym;
  rot_mtrx_RLU_G = GWinfo.symminfo.mtrx;
  for i = 1:nsym
    rot_mtrx_RLU_G{i} = double( rot_mtrx_RLU_G{i} );
  end

  % Construct rotation in Card
  rot_mtrx_Cart = cell(nsym, 1);
  rot_mtrx_RLU_R = cell(nsym, 1);
  for irot = 1:nsym
    mtrx_RLU = rot_mtrx_RLU_G{irot};
    % Since we assume all vector are row vector, operators are on r.h.s
    mtrx_Cart = a1a2a3 * mtrx_RLU * inv(a1a2a3);
    mtrx_Cart = double( mtrx_Cart );
    rot_mtrx_Cart{irot} = mtrx_Cart;
    rot_mtrx_RLU_R{irot} = (a1a2a3./fftgrid)' * mtrx_Cart * inv(a1a2a3./fftgrid)'; 
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Kpoint
  kpt_Cart = double( bz_samp.kpt );
  kpt_RLU = double( kpt_Cart / b1b2b3' );

  nibz = bz_samp.nibz;
  nsym = GWinfo.symminfo.nsym;
  nrot = GWinfo.symminfo.nrot;
  
  

  DL_vol = det( b1b2b3 );
  % 1. Construct kptbz from kpt
  nstar = int32( zeros(nibz, 1) );
  star  = int32( zeros(nibz, nsym) );
  kstar = zeros(nsym, 3);
  for ikibz = 1:nibz
    kibz_Cart = kpt_Cart(ikibz, :);
    for isym = 1:nsym
      rot_Cart = rot_mtrx_Cart{isym};
      Skibz = double( kibz_Cart * rot_Cart );
      Skibz = bz2ibz( 'Cart', Skibz, b1b2b3);
      Skibz_RLU = Skibz / b1b2b3';
      kstar(isym, :) = Skibz_RLU;
      k_found = false;
      for i_star = 1: nstar(ikibz)
        tmp = kstar(isym, :) - kstar(star(ikibz, i_star), :);
        diff = tmp - round(tmp);
        if norm( diff ) < 1e-5 * DL_vol
          k_found = true;
          break;
        end
      end
      if ~k_found
        nstar(ikibz) = nstar(ikibz) + 1;
        star(ikibz, nstar(ikibz)) = isym;
      end
    end
  end
  
  nbz = sum(nstar);
  % Construct 
  bz2ibz_ = int32( zeros(nbz, 1) );
  bz2rot = int32( zeros(nbz, 1) );
  ikbz = int32( 0 );
  for ikibz = 1:nibz
    for istar = 1:nstar(ikibz)
      ikbz = ikbz+1;
      bz2ibz_( ikbz ) = ikibz;
      bz2rot( ikbz ) = star(ikibz, istar);
    end
  end

  % Then construct kptbz
  kptbz_RLU = double( zeros(nbz, 3) );
  kptbz_Cart = double( zeros(nbz, 3) );
  for ikbz = 1:nbz
    ikibz = bz2ibz_(ikbz);
    isym = bz2rot(ikbz);
    rot_Cart = double( rot_mtrx_Cart{isym} );
    kibz_Cart = double( kpt_Cart(ikibz, :) );
    kptbz_Cart(ikbz, :) = kibz_Cart * rot_Cart;
  end
  kptbz_RLU = kptbz_Cart / b1b2b3';
  
  tmp_devel.bz2ibz = bz2ibz_;
  tmp_devel.bz2rot = bz2rot;
  tmp_devel.kptbz_Cart = kptbz_Cart;
  tmp_devel.kptbz_RLU = kptbz_RLU;
  tmp_devel.kpt_Cart = kpt_Cart;
  tmp_devel.kpt_RLU = kpt_RLU;
  tmp_devel.nibz = nibz;
  tmp_devel.nbz = nbz;
  
  % Construct qindx_S
  [qindx_S, ~] = qindx( kpt_RLU, kptbz_RLU, Ggrid_RLU, fftgrid );
  tmp_devel.qindx_S = qindx_S;


  

  % verify qindx_S
  % qindx_S = GWinfo.bz_samp.qindx_S;
  % for ikibz = 1:nibz
  %   for iqbz = 1:nbz
  %     kibz = kpt_RLU(ikibz, :);
  %     qbz = kptbz_RLU(iqbz, :);
  %     k_q = kibz - qbz;
  %     ikstar = qindx_S(ikibz, iqbz);
  %     kstar = kptbz_RLU(ikstar, :);
  %     iGo = qindx_S(ikibz, iqbz, 2);
  %     Go = Ggrid_RLU(iGo, :);
  %     % Go + kstar - (k_q)
  %   end
  % end

  % Generate inverse indices
  % rot_mtrx_*(inv_rot_index(i)) = inv(rot_mtrx_*(i))
  inv_rot_index = int32( zeros(nsym, 1)-1 );
  for irot = 1:nsym
    inv_rot_RLU_irot = inv(rot_mtrx_RLU_G{irot});
    for jrot = 1:nsym
      rot_RLU_jrot = rot_mtrx_RLU_G{jrot};
      if norm(inv_rot_RLU_irot - rot_RLU_jrot) < 1e-5
        inv_rot_index(irot) = jrot;
      end
    end
  end
  if any(inv_rot_index == -1)
    error('Error in inverse index');
  end
  inv_rot_index_RLU = inv_rot_index;

  inv_rot_index = int32( zeros(nsym, 1)-1 );
  for irot = 1:nsym
    inv_rot_Cart_irot = inv(rot_mtrx_Cart{irot});
    for jrot = 1:nsym
      rot_Cart_jrot = rot_mtrx_Cart{jrot};
      if norm(inv_rot_Cart_irot - rot_Cart_jrot) < 1e-5
        inv_rot_index(irot) = jrot;
      end
    end
  end

  %
  if any(inv_rot_index == -1)
    error('Error in inverse index');
  end
  inv_rot_index_Cart = inv_rot_index;

  if abs(inv_rot_index_RLU - inv_rot_index_Cart) > 1e-5
    error('inv_Cart .neq. inv_RLU');
  end

  tmp_devel.rot_mtrx_RLU_G = rot_mtrx_RLU_G;
  tmp_devel.rot_mtrx_RLU_R = rot_mtrx_RLU_R;
  tmp_devel.rot_mtrx_Cart = rot_mtrx_Cart;
  tmp_devel.inv_rot_index = inv_rot_index;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Generate r_rot and g_rot : nr/ng x nsym
  % Such that 
  % Rgrid(r_rot(ir, irot), :) = mtrx(irot) * Rgrid(ir, :)
  % Ggrid(g_rot(ig, irot), :) = mtrx(irot) * Ggrid(ig, :)
  R_rot = int32( zeros(nr, nsym) );
  G_rot = int32( zeros(ng, nsym) );

  for irot = 1:nsym
    mtrx_RLU = rot_mtrx_RLU_R{irot};
    if irot <= nsym / 2
      ;
    else
      mtrx_RLU = -mtrx_RLU;
    end
    rot_Rgrid_RLU = ( Rgrid_RLU * mtrx_RLU );
    for i = 1:3
      rot_Rgrid_RLU(:, i) = rot_Rgrid_RLU(:, i) - floor((rot_Rgrid_RLU(:, i)+1e-5) / fftgrid(i)) * fftgrid(i);
    end
    R_rot(:, irot) = 1 + rot_Rgrid_RLU(:, 1) + rot_Rgrid_RLU(:, 2)*fftgrid(1) ...
                   + rot_Rgrid_RLU(:, 3) * fftgrid(1) * fftgrid(2);
  end

  % G_rot
  % First, genenrate shell structure for gvec:
  normG2 = sum(Ggrid_Cart.^2, 2);
  first_index_in_each_shell = zeros(ng, 1);
  num_shell = 0;
  num_index_in_each_shell = zeros(ng, 1);
  % ig=1
  num_shell = num_shell + 1;
  num_index_in_each_shell(num_shell) = 1;
  first_index_in_each_shell(num_shell) = 1;
  %
  for ig = 2:ng
    old = normG2(ig-1);
    new = normG2(ig);
    if abs(new-old) > 1e-3
      num_shell = num_shell + 1;
      first_index_in_each_shell(num_shell) = ig;
    end
    num_index_in_each_shell(num_shell) = num_index_in_each_shell(num_shell) + 1;
  end
  % num_shell = num_shell + 1;
  % first_index_in_each_shell(num_shell) = ng;
  %   end
  % end
  num_index_in_each_shell = num_index_in_each_shell(1:num_shell);
  first_index_in_each_shell = first_index_in_each_shell(1:num_shell);
  
  for irot = 1:nsym 
    mtrx_RLU = rot_mtrx_RLU_G{irot};
    for ishell = 1:num_shell
      ig_start = first_index_in_each_shell(ishell);
      ig_len = num_index_in_each_shell(ishell);
      ig_end = ig_start + ig_len - 1;
      ig_range = ig_start:ig_end;
      rot_Ggrid_RLU = Ggrid_RLU(ig_range, :) * mtrx_RLU;
      for ig = ig_range
        G_rot(ig, irot) = find_gvec_in_glist(rot_Ggrid_RLU(ig-ig_start+1, :), ...
                          Ggrid_RLU(ig_range, :), fftgrid, 1e-5) + ig_start-1;
      end
    end
  end
  
  % compare with GWinfo.mapping.g_rot
  for irot = 1:nsym
    if ~isempty(GWinfo.mapping.g_rot{irot})
      perm1 = double( round( GWinfo.mapping.g_rot{irot} ) );
      perm2 = double( round( G_rot(:, irot) ) );
      fprintf("irot = %3d, diff = %.3e.\n", irot, norm(perm1 - perm2));
    end
  end

  tmp_devel.R_rot = R_rot;
  tmp_devel.G_rot = G_rot;
  tmp_devel.num_shell = num_shell;
  tmp_devel.first_index_in_each_shell = first_index_in_each_shell;
  tmp_devel.num_index_in_each_shell = num_index_in_each_shell;
  % find Gomapping, then compare with setGomap --> GWinfo.Gomapping
  % Assume well-done

  % Construct 
  % Construct g_table : G - G_0 --> G
  % max_Go = max( GWinfo.b)
  % for ig2 = 1:max()

  % Construct the full BZ sampling
  % for ik_ibz = 1, nibz
  %   for i = 1:nsym
  %     v =  





  
  
  GWinfo.tmp_devel = tmp_devel;
  
end

function [qindx_S, qindx_X] = qindx(kpt_ibz, kpt_bz, Ggrid, fftgrid)
  % qindx_S(ikibz, iqbz, 1) = ikpbz
  % qindx_S(ikibz, iqbz, 2) = iGo
  % kp + Go = k - q
  % 
  nibz = size(kpt_ibz, 1);
  nbz = size(kpt_bz, 1);
  
  %
  qindx_S = int32( zeros(nibz, nbz, 2) );
  qindx_X = int32( zeros(nbz, nbz, 2) );
  %
  for ikibz = 1:nibz
    kibz = kpt_ibz(ikibz, :);
    for iqbz = 1:nbz
      qbz = kpt_bz(iqbz, :);
      k_q = kibz - qbz;
      for ikpbz = 1:nbz
        kpbz = kpt_bz(ikpbz, :);
        kpt_diff = k_q - kpbz;
        if norm( kpt_diff - round(kpt_diff) ) < 1e-5
          break;
        end
        if ikpbz == nbz
          error('No matching kpbz found');
        end
      end
      qindx_S(ikibz, iqbz, 1) = ikpbz;
      
      % 
      G0 = kibz - qbz - kpbz;
      G0 = round(G0);
      indG0 = find_gvec_in_glist( G0, Ggrid, fftgrid );
      qindx_S(ikibz, iqbz, 2) = indG0;
    end
  end

end