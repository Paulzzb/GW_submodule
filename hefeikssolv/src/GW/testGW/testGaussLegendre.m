addpath(genpath("../../../"));
f = @(x) exp(sqrt(-1) * x);
nmax = 20;

n = nmax+5;
% GaussLegendre_;
% nodes = theta;
% weights = s;

reslist = zeros(nmax, 1)

for n = 1:nmax
% Calculate 
%   \int_0^\infty f(x) dx 
% = \int_{-1}^1 f(x(z)) dx(z)
% = \int_{-1}^1 f(x(z)) * 2 / (1-z)^2 dz
% \approx \sum_{k=1}^n c_k * F(z_k),
% where 
% {c_k}, {z_k} are the Gauss quadrature weights and quadrature nodes
%              on [-1, 1] with weight function w == 1.
% x(z) = 2/(1-z) - 1;
% F(z) = f(x(z)) * 2 / (1-z)^2

  [s, theta] = GaussLegendre(n);
  % Now, c_k = s(k) and z_k = theta(k)
  znodes = theta;
  xz = @(z) (2 ./ (1-z) -1);
  xnodes = xz(znodes);
  fxz = f(sqrt(-1)*xnodes);
  
  
  % Calculate F(z)
  Fz =  fxz * 2 ./ (1-znodes).^2;
  
  intGauss = sum(Fz .* s);
  fprintf("nodes = %d, intGauss = %12.6f, diff = %12.6f.\n", ...
         n, intGauss, intGauss-1);

  % Calculate reslist
  % reslist(n) = \int_k=1 ^ n \Pi_{k=1}^n (x-x_k)^2 dx
  poly = 1; polytmp = 0;
  for i = 1:n
    polytmp = conv([1, -znodes(i)], [1, -znodes(i)]);
    poly = conv(poly, polytmp);
  end 
  reslist(n) = polyint(poly);
end

% function [, theta]
function out = polyint(p1, p2)
  if nargin < 2
    p2 = 1;
  end
  p1 = reshape(p1, 1, []);
  p2 = reshape(p2, 1, []);
  p = conv(p1, p2);
  degp = length(p);
%   out = (degp:-1:1)';
%   out = sum(p ./ out);
  out = zeros(1, degp);
  for i = degp:-1:1
    deg = degp - i;
    if (mod(deg, 2) ~= 0)
      out(i) = 0;
    else
      out(i) = 2 ./ (deg+1);
    end
  end
  out = sum(p .* out);
end

