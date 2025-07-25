function [bz_samp] = fullbz(bz_samp, syms, gvec)
TOL = 1e-9;
nibz = bz_samp.nibz;
kbz_ca = zeros(syms.nrot,3);
kibz_ca = round(bz_samp.kpt / bz_samp.bmatrix, 6); % reduce accuracy to prevent truncation errors
nbz = 0;
for ir=1:nibz
  for it=1:syms.ntran
    tmpf = kibz_ca(ir,:) * syms.mtrx{it,1};
    [tmpf, gpt] = krange(tmpf,TOL);
    found = 0;
    for ifull=1:nbz
      if all(abs(tmpf-kbz_ca(ifull,:)) < TOL)
        found = 1;
        break
      end
    end
    if (found == 1)
      continue
    else
      nbz = nbz+1;
      kbz_ca(nbz,:)=tmpf;
      ind_rotation(nbz)=it;
      ind_kbz(nbz)=ir;
      gptlist(nbz, 1:3)=gpt;
    end
  end
end

kbz_ca = kbz_ca(1:nbz, :);
bz_samp.kptbz = kbz_ca(1:nbz, :) * bz_samp.bmatrix;
bz_samp.nbz = nbz;
bz_samp.kbz2kibz_ind_rotation=ind_rotation;
bz_samp.kbz2kibz_ind_kbz=ind_kbz;



% Generate nGo based on kbz2kibz_ind_G0
if (gvec.ng >= 27)
  maxGoset = gvec.components(1:27, :);
else
  warning('gvec.ng < 27, might need bigger ecut for ');
  maxGoset = gvec.components;
end
iGolist = zeros(nbz, 1);

for ibz = 1:nbz
  indGo = find_vec_in_list(gptlist(ibz, :), maxGoset, TOL);
  iGolist(ibz) = indGo;
end
bz_samp.nGo = max(iGolist);
bz_samp.iGolist = iGolist;




if (true)
  for ii=1:nbz
    for jj=1:ii-1
      tmpf=abs(kbz_ca(ii,:)-kbz_ca(jj,:));
      tmpf=tmpf-floor(tmpf);
      tmpf(tmpf >= 0.5) = 1 - tmpf(tmpf >= 0.5);
      if sum(abs(tmpf)) <= 1e-9
        error('equivalent points found in the full BZ, equiv kpts %d and %d with diff %d', ii, jj, kptbz(ii,:)-kptbz(jj,:));
      end
    end
  end
end


% Second, calculate qindx_* in bz_samp
qindx_S = zeros(nibz, nbz, 2);
for ik = 1:bz_samp.nibz
  for iqbz = 1:bz_samp.nbz
    % ikbz_ibz = bz_samp.kbz2kibz_ind_kbz(ikbz);

    % Find corresponding okbz = ikbz-iqbz and 
    qpt = kbz_ca(iqbz,:);
    kptbz = kibz_ca(ik,:);
    k_q = kptbz - qpt;
    [k_qbz, g0] = krange(k_q, TOL);
    ind_k = find_vec_in_list(k_qbz, kbz_ca);
    ind_g0 = find_vec_in_list(g0, maxGoset);
    % ind_g0 = find_vec_in_list(g0, bz_samp.Go_list);
    if ind_k == -1
      error('k_qbz not found in kptbz_ca');
    end
    if ind_g0 == -1
      error('g0 not found in Go_list');
    end
    qindx_S(ik, iqbz, 1) = ind_k;
    qindx_S(ik, iqbz, 2) = ind_g0;
    % bz_samp.qindx_X(iq, ikbz) = g0;
  end
end

qindx_X = zeros(nibz, nbz, 2);
for iq = 1:bz_samp.nibz
  for ikbz = 1:bz_samp.nbz
    % ikbz_ibz = bz_samp.kbz2kibz_ind_kbz(ikbz);

    % Find corresponding okbz = ikbz-iqbz and 
    qpt = kibz_ca(iq,:);
    kptbz = kbz_ca(ikbz,:);
    k_q = kptbz - qpt;
    [k_qbz, g0] = krange(k_q, TOL);
    ind_k = find_vec_in_list(k_qbz, kbz_ca);
    ind_g0 = find_vec_in_list(g0, maxGoset);
    % ind_g0 = find_vec_in_list(g0, bz_samp.Go_list);
    if ind_k == -1
      error('k_qbz not found in kptbz_ca');
    end
    if ind_g0 == -1
      error('g0 not found in Go_list');
    end
    qindx_X(iq, ikbz, 1) = ind_k;
    qindx_X(iq, ikbz, 2) = ind_g0;
    % bz_samp.qindx_X(iq, ikbz) = g0;
  end
end
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




