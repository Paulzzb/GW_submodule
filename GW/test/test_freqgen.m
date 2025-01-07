% Using g(z) = exp(iz) to calculate 
%    \int_0^{i+\infty} f(z)g(z) dz
% where f(z;)|_{real(z) = 0} = a / (a^2+|z|^2).
addpath(genpath("../../../"));

ii = sqrt(-1);
a = 3;
f = @(z, a) a ./ (a.^2 + abs(z).^2);
g = @(z) exp(ii*z).*(1+abs(z).^2+abs(z).^4);
g_ = @(z) exp(-z).*(1+abs(z).^2+abs(z).^4);
xz = @(z) (2 ./ (1-z) -1);

n = 20;
[s, theta] = GaussLegendre(n);
dx = 2 ./ (1-theta).^2;
s = s .* dx;
theta = ii * xz(theta);

result1 = sum(s .* f(theta, a) .* g(theta));

N = 10;
h = 1e-3;
nodes = ii * ((0:h:(N-h)));
nnodes = length(nodes);
weights = h * ones(1, nnodes);

result2 = sum(weights .* f(nodes, a) .* g(nodes));

h = @(x) f(x, a).*g_(x);
result3 = integral(h, 0, N);

result3 - result1
result3 - result2

% Reconstruct how BGW calculate
n4 = 100;
h = 1/n4;
zfreq = 0:h:1-h;
zfreq = zfreq ./ (1-zfreq);
zfreqmid = zeros(1, n4); zfreqmid(1) = zfreq(1);
for i = 2:n4
  zfreqmid(i) = (zfreq(i-1)+zfreq(i))/2;
end
zfreq = zfreq * ii;
zfreqmid = zfreqmid * ii;

F = @(x, a) atan(abs(x) / a);
result4 = 1/2 * (F(zfreqmid(2), a) - F(0, a)) * g(zfreq(1));
for i = 2:n4-1
  result4 = result4 + (F(zfreqmid(i+1), a) - F(zfreqmid(i), a)) * g(zfreq(i));
end
result4 = result4 + 0.5 * (F(zfreq(n4), a) - F(zfreq(n4), a)) * g(zfreq(i));

result2 - result4
