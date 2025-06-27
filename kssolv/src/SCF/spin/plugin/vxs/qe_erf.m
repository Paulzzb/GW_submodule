function erf =  qe_erf(x)
% x must be real!
p1 =  [2.426679552305318E2 2.197926161829415E1...
    6.996383488619136 -3.560984370181538E-2];
q1 =  [2.150588758698612E2 9.116490540451490E1...
    1.508279763040779E1 1.000000000000000];
ax = abs(x);

erf = zeros(size(x));

idxl = ax > 6.0;
erf(idxl) = sign(x(idxl));

idxs = ax <= 0.47; x2 = x(idxs).^2;
erf(idxs) = x(idxs) .* (p1 (1) + x2 .* (p1 (2) + x2 .* (p1 (3) + x2 * p1 (4) ) ) ) ...
    ./ (q1 (1) + x2 .* (q1 (2) + x2 .* (q1 (3) + x2 * q1 (4) ) ) );

idxm = ~idxl & ~idxs;
if nnz(idxm)
    erf(idxm) = 1 - qe_erfc(x(idxm));
end

end