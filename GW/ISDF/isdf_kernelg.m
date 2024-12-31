function zetag = isdf_kernel(Phi, Psi, ind_mu, gvec, vol)

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

if m1 ~= m2,
    error('Wrong inputs: row dimensions of Phi and Psi do not match!\n');
end
m = m1;
if length(ind_mu) > m,
    error('Wrong inputs: ind_mu is too long!\n');
end

rk=length(ind_mu);
phi=Phi(ind_mu,:);
psi=Psi(ind_mu,:);
C2=(phi*phi').*(psi*psi');
% L = chol(C2, 'lower');
% invL = inv(L);


% Slicing C1 if m is too big.
rank_mu = length(ind_mu);
% if m*rank_mu < 1e+9
if 0
  C1=(Phi*phi').*(Psi*psi');
  zetar=C1/C2;
  clear C1 C2;
  parfor i = 1:rk
    fftbox1 = reshape(zetar(:, i), gvec.fftgrid);
    fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * vol;
    zetag(:, i) = get_from_fftbox(gvec.idxnz, fftbox1, gvec.fftgrid);
  end
else
  % C1g = zeros(ng, rk);
  % imax = ceil(rank_mu / step);
  % for i = 1:imax
  %   if i == 2
  %     startfirstiter = tic;
  %   end
  %   if i < imax
  %     irange = (i-1)*step+1:i*step;
  %   else
  %     irange = (i-1)*step+1:imax;
  %   end
  %   tmpr = (Phi * phi(irange, :)') .* (Psi * psi(irange, :)');
  %   for j = 0:length(irange)-1
  %     fftbox1 = reshape(tmpr(:, j+1), gvec.fftgrid);
  %     fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * vol;
  %     C1g(:, i+j) = get_from_fftbox(gvec.idxnz, fftbox1, gvec.fftgrid);
  %   end
  %   if i == 5*step+1
  %     timefirst = toc(startfirstiter);
  %     fprintf("%20s = %6.2f", 'Time estimate for isdf_kernel.\n', timefirst / (3*step) * rk);
  %   end
  % end
  % zetag = C1g / C2;
  C1g = zeros(gvec.ng, rk);
  for i = 1:step:rk
    if i == 2*step+1
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
    if i == 5*step+1
      timefirst = toc(startfirstiter);
      fprintf("%20s = %6.2f", 'Time estimate for isdf_kernel.\n', timefirst / (3*step) * rk);
    end
  end
  zetag = C1g / C2;
end

return;
