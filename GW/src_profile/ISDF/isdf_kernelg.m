function zetag = isdf_kernelg(Phi, Psi, ind_mu, gvec, vol)

% This file will be final version for the full-frequency GW calculation.
% 2024/11/07           Zhengbang Zhou     Fudan University.

% Description
%     Using indices from ISDF main code to generate helper functions
%     in G-space directly. More specifically, we are calculating the
%     following problem:
%      zetag = iFFT(  (Phi*phi').*(Psi*psi') / (phi*phi').*(psi*psi')  )
%     where phi = Phi(ind_mu, :), psi = Psi(ind_mu, :).
% 


% Parameters
%     Input:
%         Phi, Psi: original functions.
%         ind_mu: indices from ISDF main code.
%         gvec: class Ggrid, which is used to perform FFT.
%         vol: volumn of the system in real space, which is used to perform FFT.
%     Output:
%         zetag_mu: helper functions in G-space.

% Structure 
%     The code is organized as followed:
%     1. Calculate C2 = (phi*phi').*(psi*psi'); 
%     2. Calculate C1 = (Phi*phi').*(Psi*psi');
%     3. Do iFFT to C1, C1g = iFFT(C1);
%     4. Calculate helper functions to zetag = C1g / C2;


[m1, n1] = size(Phi);
[m2, n2] = size(Psi);
ng = gvec.ng;
step = 25;

if m1 ~= m2
  msg = 'Wrong inputs: row dimensions of Phi and Psi do not match!\n';
  GWerror(msg);
end
m = m1;
if length(ind_mu) > m
  msg = 'Wrong inputs: ind_mu is too long!\n';
  GWerror(msg);
end

rk=length(ind_mu);
phi=Phi(ind_mu,:);
psi=Psi(ind_mu,:);
C2=(phi*phi').*(psi*psi');


fprintf('\n[ISDF] Generating helper functions in G-space...\n');
rank_mu = length(ind_mu);
total_iter = ceil(rk / step);
est_total_time = -1;


C1g = zeros(gvec.ng, rk);



for i = 1:step:rk
  iter_idx = ceil(i / step);

  if iter_idx == 2
    startfirstiter = tic;
  end

  if i+step < rk
    irange = i:i+step-1;
  else
    irange = i:rk;
  end

  tmpr = (Phi * phi(irange, :)') .* (Psi * psi(irange, :)');
  for j = 0:length(irange)-1
    fftbox1 = reshape(tmpr(:, j+1), gvec.fftgrid);
    fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * vol;
    C1g(:, i+j) = get_from_fftbox(gvec.idxnz, fftbox1, gvec.fftgrid);
  end

  if iter_idx == 6 
    time_per_step = toc(startfirstiter) / 4;
    est_total_time = time_per_step * total_iter;
    fprintf('[ISDF helper] Estimated total time: %.1f seconds\n', est_total_time);
  end

  % Progress bar
  if est_total_time > 0
    elapsed_time = (iter_idx - 1) * time_per_step;
    fprintf('[ISDF helper] Progress: %3d%% | Elapsed: %.1fs / Estimated: %.1fs\n', ...
        round(100 * (iter_idx-1) / total_iter), elapsed_time, est_total_time);
  else
    fprintf('[ISDF helper] Progress: %3d%%\n', round(100 * (iter_idx-1) / total_iter));
  end
end

zetag = C1g / C2;

return;
