function zetag = isdf_kernel(Phi, Psi, ind_mu, gvec, vol)

% Using indices from ISDF main code to generate helper functions
% in G-space.
% This code is designed to reduce the memory cost when generating
% helper functions in R-space, which dominate the memory cost when
% the system is relatively big.
% Input:
%   Phi, Psi: original functions.
%   ind_mu: indices from ISDF main code.
%   gvec: class Ggrid, which is used to perform FFT.
%   vol: volumn of the system in real space, which is used to perform FFT.
% Output:
%   zetag_mu: helper functions in G-space.

[m1, n1] = size(Phi);
[m2, n2] = size(Psi);

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
  for i = 1:rk
    phir = (Phi * phi(i, :)') .* (Psi * psi(i, :)');
    fftbox1 = reshape(phir, gvec.fftgrid);
    fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * vol;
    C1g(:, i) = get_from_fftbox(gvec.idxnz, fftbox1, gvec.fftgrid);
  end
  zetag = C1g / C2;
end

return;
