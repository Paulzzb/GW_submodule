function [bz_samp] = fullbz(bz_samp, syms, gvec)

TOL = 1e-4;
nibz = bz_samp.nibz;
kibz_RLU = round(bz_samp.kpt / bz_samp.bmatrix, 6); % reduce accuracy to prevent trunRLUtion errors
nbz = 0;
fftgrid = gvec.fftgrid; 
bmatrix = bz_samp.bmatrix;
nrot = syms.nrot;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Kpoint
  b1b2b3 = bmatrix';
  kpt_Cart = single( bz_samp.kpt );
  kpt_RLU = single( kpt_Cart / b1b2b3' );

  nibz = bz_samp.nibz;
  nrot =  syms.nrot;
  ntran = syms.nrot;
  rot_mtrx_RLU = syms.mtrx;
  for isym = 1:nrot
    rot_mtrx_Cart{isym} = single( inv(b1b2b3)' * syms.mtrx{isym} * b1b2b3' );
  end
  
  
  

  DL_vol = det( b1b2b3 );
  % 1. Construct kptbz from kpt
  nstar = int32( zeros(nibz, 1) );
  star  = int32( zeros(nibz, nrot) );
  kstar = zeros(nrot, 3);
  for ikibz = 1:nibz
    kibz_Cart = kpt_Cart(ikibz, :);
    for isym = 1:nrot
      rot_Cart = rot_mtrx_Cart{isym};
      Skibz = single( kibz_Cart * rot_Cart );
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
  kptbz_RLU = single( zeros(nbz, 3) );
  kptbz_Cart = single( zeros(nbz, 3) );
  for ikbz = 1:nbz
    ikibz = bz2ibz_(ikbz);
    isym = bz2rot(ikbz);
    rot_Cart = single( rot_mtrx_Cart{isym} );
    kibz_Cart = single( kpt_Cart(ikibz, :) );
    kptbz_Cart(ikbz, :) = kibz_Cart * rot_Cart;
  end
  kptbz_RLU = kptbz_Cart / b1b2b3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kbz_RLU = kptbz_RLU;
ind_rotation = bz2rot;
ind_kbz = bz2ibz_;
gptlist = zeros(nbz, 3);




% for ir=1:nibz
%   for it=1:syms.ntran
%     Skibz = kibz_RLU(ir,:) * syms.mtrx{it,1};
%     % [tmpf, gpt] = krange(tmpf,TOL);
%     tmpf = bz2ibz( 'rlu', Skibz, bz_samp.bmatrix');
%     gpt = round(tmpf - Skibz);
%     found = 0;
%     for ifull=1:nbz
%       diff = tmpf - kbz_RLU(ifull, :);
%       if all(abs( tmpf-kbz_RLU(ifull,:) ) < TOL)
%         found = 1;
%         break
%       end
%     end
%     if (found == 1)
%       continue
%     else
%       nbz = nbz+1;
%       kbz_RLU(nbz,:)=tmpf;
%       ind_rotation(nbz)=it;
%       ind_kbz(nbz)=ir;
%       gptlist(nbz, 1:3)=gpt;
%     end
%   end
% end

% We need to remove redundant points from the full BZ.
% In time_symm, then -I_3 is in the group, so it is possible, that k and -k are both in the full BZ.
% However, this is not allowed, thus we need to 
% if (true)
%   for ii=1:nbz
%     for jj=1:ii-1
%       tmpf=abs(kbz_RLU(ii,:)-kbz_RLU(jj,:));
%       tmpf=tmpf-floor(tmpf);
%       tmpf(tmpf >= 0.5) = 1 - tmpf(tmpf >= 0.5);
%       if sum(abs(tmpf)) <= 1e-9
%         error('equivalent points found in the full BZ, equiv kpts %d and %d with diff %d', ii, jj, kbz_RLU(ii,:)-kbz_RLU(jj,:));
%       end
%     end
%   end
% end



kbz_RLU = kbz_RLU(1:nbz, :);
bz_samp.kptbz = kbz_RLU(1:nbz, :) * bz_samp.bmatrix;
bz_samp.nbz = nbz;
bz_samp.kbz2kibz_ind_rotation=ind_rotation;
bz_samp.kbz2kibz_ind_kbz=ind_kbz;



% Generate nGo based on kbz2kibz_ind_G0
% if (gvec.ng >= 27)
%   maxGoset = gvec.components(1:27, :);
% else
%   warning('gvec.ng < 27, might need bigger ecut for ');
  maxGoset = gvec.components;
% end

iGolist = zeros(nbz, 1);

for ibz = 1:nbz
  indGo = find_gvec_in_glist(gptlist(ibz, :), maxGoset, fftgrid, TOL);
  iGolist(ibz) = indGo;
end
bz_samp.nGo = max(iGolist);
bz_samp.iGolist = iGolist;







% Second, RLUlculate qindx_* in bz_samp
qindx_S = zeros(nibz, nbz, 2);
for ik = 1:bz_samp.nibz
  for iqbz = 1:bz_samp.nbz
    % ikbz_ibz = bz_samp.kbz2kibz_ind_kbz(ikbz);

    % Find corresponding okbz = ikbz-iqbz and 
    iq_ibz = ind_kbz(iqbz);
    iq_s = ind_rotation(iqbz);
    qpt = kbz_RLU(iqbz,:);
    qpt_ibz = kibz_RLU(iq_ibz, :);
    Sqpt = syms.mtrx{iq_s, 1};
    %
    kptbz = kibz_RLU(ik,:);
    Sqibz = qpt_ibz * Sqpt;
    k_Sqibz = kptbz - Sqibz;
    kp_bz = bz2ibz( 'rlu', k_Sqibz, bz_samp.bmatrix');
    % [kp_bz, ~] = krange(k_Sqibz, TOL);
    % g0 = -g0;
    ind_kp = find_kvec_in_klist(kp_bz, kbz_RLU);
    % ind_g0 = find_gvec_in_glist(g0, maxGoset, fftgrid);
    % ind_g0 = find_vec_in_list(g0, bz_samp.Go_list);
    % k_Sq_kp = k_Sqibz - kp_bz;
    k_Sq_kp = kptbz - qpt - kp_bz;
    % k_Sq_kp = -k_Sq_kp;
    disp([kptbz; Sqibz; qpt; kp_bz])
    ind_g0 = find_gvec_in_glist(k_Sq_kp, maxGoset, fftgrid);
    if ind_kp == -1
      error('k_qbz not found in kptbz_RLU');
    end
    if ind_g0 == -1
      error('g0 not found in Go_list');
    end
    qindx_S(ik, iqbz, 1) = ind_kp;
    qindx_S(ik, iqbz, 2) = ind_g0;
    % bz_samp.qindx_X(iq, ikbz) = g0;
  end
end

qindx_X = zeros(nibz, nbz, 2);
% for ik_ibz = 1:bz_samp.nibz
%   for iq_bz = 1:bz_samp.nbz
%     % ikbz_ibz = bz_samp.kbz2kibz_ind_kbz(ikbz);
% 
%     % Find corresponding okbz = ikbz-iqbz and 
%     kpt_ibz = kibz_RLU(ik_ibz,:);
%     qpt_bz = kbz_RLU(iq_bz,:);
%     k_q = kpt_ibz - qpt_bz;
%     [k_qbz, g0] = krange(k_q, TOL);
%     g0 = -g0;
%     ind_k = find_kvec_in_klist(k_qbz, kbz_RLU);
%     ind_g0 = find_gvec_in_glist(g0, maxGoset, fftgrid);
%     % ind_g0 = find_vec_in_list(g0, bz_samp.Go_list);
%     if ind_k == -1
%       error('k_qbz not found in kptbz_RLU');
%     end
%     if ind_g0 == -1
%       error('g0 not found in Go_list');
%     end
%     qindx_X(ik_ibz, iq_bz, 1) = ind_k;
%     qindx_X(ik_ibz, iq_bz, 2) = ind_g0;
%     % bz_samp.qindx_X(iq, ikbz) = g0;
%   end
% end
bz_samp.qindx_S = qindx_S;
bz_samp.qindx_X = qindx_X;

end % EOF

% function [ind_k] = find_vec_in_list(kpt, kpt_set, TOL)
%   % Return the index of kpt in kpt_set
%   %       -1 if not found
%   % Default TOL = 1e-9
%   if nargin < 3
%     TOL = 1e-9;
%   end

%   ind_k = -1;
%   nbz = length(kpt_set(:)) / 3;

%   for ii=1:nbz
%     tmpf=abs(kpt-kpt_set(ii,:));
%     if sum(abs(tmpf)) <= TOL 
%       ind_k = ii;
%       return;
%     end
%   end

%   warning("find_vec_in_list: kpt not found in kpt_set");
% end % EOF 




