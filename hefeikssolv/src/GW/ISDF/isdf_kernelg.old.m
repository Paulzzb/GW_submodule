function zetag = isdf_kernelg(Phi, Psi, ind_mu, gvec, vol)
% This file will be final version for the full-frequency GW calculation.

% Description
%     Using indices from ISDF main code to generate helper functions
%     in G-space directly.
%     Compare to normal process, where first get the helper functions in
%     R-space, then transform it to G-space, this reduce memory cost when
%     the system is relatively big.

% Parameters
%     Input:
%         Phi, Psi: original functions.
%         ind_mu: indices from ISDF main code.
%         gvec: class Ggrid, which is used to perform FFT.
%         vol: volumn of the system in real space, which is used to perform FFT.
%     Output:
%         zetag_mu: helper functions in G-space.

[m1, n1] = size(Phi);
[m2, n2] = size(Psi);

if m1 ~= m2,
    error('Wrong inputs: row dimensions of Phi and Psi do not match!\n');
end
m = m1;
if length(ind_mu) > m,
    error('Wrong inputs: ind_mu is too long!\n');
end
Phigpu = gpuArray(Phi);
Psigpu = gpuArray(Psi);

rk=length(ind_mu);
nsteprows = 25;
stepcols = 10000;
phi=Phi(ind_mu,:);
psi=Psi(ind_mu,:);
phigpu = gpuArray(phi);
psigpu = gpuArray(psi);
C2= (phi*phi').*(psi*psi');
L = chol(C2, 'lower');
invL = inv(L);

% Slicing C1 if m is too big.
rank_mu = length(ind_mu);
% if m*rank_mu < 1e+9
if 0
  C1=(Phi*phi').*(Psi*psi');
  zetar=C1/C2;
  clear C1 C2;
  for i = 1:rk
    fftbox1 = reshape(zetar(:, i), gvec.fftgrid);
    fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * vol;
    zetag(:, i) = get_from_fftbox(gvec.idxnz, fftbox1, gvec.fftgrid);
  end
else
  C1g = zeros(gvec.ng, rk);
  for itmp = 1:rk/step
    i = 1 + step*(itmp-1)
    if i == 2*step+1
      startfirstiter = tic;
    end
    if i+step < rk
      irange = i:i+step-1;
    else
      irange = i:rk;
    end
    tmprgpu = (Phigpu * phigpu(irange, :)') .* (Psigpu * psigpu(irange, :)');
    for j = 0:length(irange)-1
      fftboxgpu = reshape(tmprgpu(:, j+1), gvec.fftgrid);
      % fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * vol;
      fftbox1gpu = ifftn(fftboxgpu, gvec.fftgrid);
      % tmpgpu = get_from_fftbox(gvec.idxnz, fftbox1gpu, gvec.fftgrid);
      tmpgpu = fftbox1gpu(gvec.idxnz);
      C1g(:, i+j) = gather(tmpgpu);
    end
    if i == 5*step+1
      timefirst = toc(startfirstiter);
      fprintf("%20s = %6.2f", 'Time estimate for isdf_kernel.\n', timefirst / (3*step) * rk);
    end
  end
  clear tmpgpu, fftbox1gpu, fftboxgpu;
  ng = size(C1g, 1);
  for itmp = 1:ng/step
    i = 1+(itmp-1)*step;
    if i+step < ng
      irange = i:i+step-1;
    else
      irange = i:ng;
    end
    C1ggpu = gpuArray(C1g(irange, :));
    tmpgpu = C1ggpu / C2gpu;
    zetag(irange, :) = gather(tmpgpu);
  end
end

return;
