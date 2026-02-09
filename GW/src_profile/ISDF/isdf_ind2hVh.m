% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/02/02 ZZ
function hVh = isdf_ind2hVh(Phi, Psi, ind_mu, Dcoul, gvec, vol)
% Formly, we have helper = (\F*M C')*(CC')^{-1}, hVh = <h | Dcoul | h>
% So, we calculate it as
% 1. Column slicing to calculate TMP = \F*(M C') 
% 2. Calculate TMP2 = TMP'*Dcoul*TMP;
% 3. Calculate C2 = CC', TMP4 = pseudo_inv(CC')? (or direct inverse)
% 4. Calculate hVh = TMP4' * TMP2 * TMP4

cleanup = QPlog_push('ISDF-ind2hVh');
msg = sprintf('Log to be done');
QPlog(msg, 2);

[m1, n1] = size(Phi);
[m2, n2] = size(Psi);
ng = gvec.ng;
Nblock = 64;
% Nblock = 16*gcp('nocreate').NumWorkers;

if m1 ~= m2
  msg = 'Wrong inputs: row dimensions of Phi and Psi do not match!';
  QPerror(msg);
end
m = m1;
if length(ind_mu) > m
  msg = 'Wrong inputs: ind_mu is too long!';
  QPerror(msg);
end

rk=length(ind_mu);
phi=Phi(ind_mu,:);
psi=Psi(ind_mu,:);
fftgrid = gvec.fftgrid;
idxnz = gvec.idxnz;

% -------------------------------------------------------------------
% 3. Calculate C2 = CC', TMP4 = pseudo_inv(CC')? (or direct inverse)
C2=(phi*phi').*(psi*psi');
if condest(C2) < 1e+12
  C2inv = inv(C2);
else
  C2inv = pseudoinv(C2);
end
clear C2;

% -------------------------------------------------------------------
% 1. Column slicing to calculate TMP = \F*(M C') 
total_iter = ceil(rk / Nblock);
est_total_time = -1;

C1g = zeros(gvec.ng, rk);
for i = 1:Nblock:rk
  iter_idx = ceil(i / Nblock);

  if iter_idx == 2
    startfirstiter = tic;
  end

  if i+Nblock < rk
    irange = i:i+Nblock-1;
  else
    irange = i:rk;
  end

  tmpr = (Phi * phi(irange, :)') .* (Psi * psi(irange, :)');
  for j = 0:length(irange)-1
  % parfor j = 0:length(irange)-1
    fftbox1 = reshape(tmpr(:, j+1), fftgrid);
    fftbox1 = do_FFT(fftbox1, fftgrid, 1) * vol;
    C1g(:, i+j) = get_from_fftbox(idxnz, fftbox1, fftgrid);
  end

  if iter_idx == 6 
    time_per_Nblock = toc(startfirstiter) / 4;
    est_total_time = time_per_Nblock * total_iter;
    msg = sprintf('Estimated total time: %.1f seconds', est_total_time);
    QPlog(msg, 2);
  end

  % Progress bar
  if est_total_time > 0
    elapsed_time = (iter_idx - 1) * time_per_Nblock;
    msg = sprintf('Progress: %3d%% | Elapsed: %.1fs / Estimated: %.1fs', ...
        round(100 * (iter_idx-1) / total_iter), elapsed_time, est_total_time);
    QPlog(msg, 2);
  else
    msg = sprintf('Progress: %3d%%', round(100 * (iter_idx-1) / total_iter));
    QPlog(msg, 2);
  end
end
% -------------------------------------------------------------------
% 2. Calculate TMP2 = TMP'*Dcoul*TMP;
C1VC1 = zeros(rk, rk);
for i = 1:Nblock:rk
  iter_idx = ceil(i / Nblock);
  if iter_idx == 2
    startfirstiteri = tic;
  end
  if i+Nblock < rk
    irange = i:i+Nblock-1;
  else
    irange = i:rk;
  end
  %
  for j = 1:Nblock:rk
    iter_jdx = ceil(j / Nblock);
    if iter_jdx == 2
      startfirstiterj = tic;
    end

    if j+Nblock < rk
      jrange = j:j+Nblock-1;
    else
      jrange = j:rk;
    end
    %
    C1VC1(irange, jrange) = C1g(:, irange)' * Dcoul * C1g(:, jrange) / vol;
  end
end

% -------------------------------------------------------------------
% 4. Calculate hVh = TMP4' * TMP2 * TMP4
hVh = zeros(rk, rk);
hVh = C2inv * C1VC1 * C2inv;

end % function