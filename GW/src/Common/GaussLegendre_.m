function [s, theta] = GaussLegendre(n, tol)
%     [s, theta] = GaussLegendre(n, tol)
% Generate the Gauss quadrature rule on [-1, 1]
% Description:
% Gauss--Legendre quadrature rule generate a set of quadrature weights
% {c_k} and quadrature nodes {x_k}, and the integral can be approximated
% with 
% \int_{-1}^1 f(x) dx = \sum_{k=1}^n c_k f(x_k) + r_n(f),
% where the residual r_n(f) is bound by 
% |r_n(f)| <= f^(2n)(x)/(2*n!) \int_{-1}^1 \multi_{k=1}^n (x-x_k)^2 dx. 
% 
% Since we use the Gauss--Legendre quadrature rule to deal with
      % \int_{0}^{i*\infty} G(w+w')W(w') dw',
% and G(w+w')W(w') acts like 1 / (w'-a) on imaginary axis.
% So in this particular situation, we expect that 
      % |r_n(f)| = O( \int_{-1}^1 \multi_{k=1}^n (x-x_k)^2 dx). 
% Based on that estimation, we generate a "self-adaption method" with
% respect to the required tolerate 'tol'.
% Input:
%     n: number of quadrature nodes.
%     tol (optional): required tolerate,
%                     if given, then 'n' is rewritten by a new one.
% Output:
%     s: weights
%     theta: nodes

TOL_ZERO = 1e-12;

if nargin == 2
  if tol < 1e-10
    fprintf("Tolerate for integral is setted too small: %8.3e.\n", tol);
  end
  % Here, the reslist is used to decided which n to choose from
  reslist = [
       6.666e-01
       1.777e-01
       4.571e-02
       1.160e-02
       2.931e-03
       7.380e-04
       1.854e-04
       4.654e-05
       1.167e-05
       2.925e-06
       7.329e-07
       1.835e-07
       4.595e-08
       1.150e-08
       2.879e-09
       7.294e-10
       1.826e-10
       5.367e-11
       4.169e-11
       2.183e-11];
  
  for n = 1:length(reslist)
    if reslist(n) < tol
      break;
    end
  end
end

if n > 20
  warning(["Gauss quadrature rule is not stable with n = %3d.\n",...
           "Reset to be %3d for safety.\n"], n, 15);
  n = 15;
end


% Generate the Gauss Quadrature rule on [0, 1] with n nodes.
pi = zeros(1, 2*n);
pi_1 = 1; % pi_1(end) = 1; % p_0 = 1
pi_2 = 0; % p_{-1} = 0;
x = [1, 0]; % represent polynormial p(x) = x;
T = zeros(n, n);
beta = 0;
alpha = 0;


% Generate T using three terms iteration
%   beta_i p_i = (x - alpha_i)p_{i-1} - beta_{i-1}p_{i-2} = tilde_pi
% where
%   beta_i  : normalizing factor, beta_i = sqrt(<tildepi, tildepi>)
%   alpha_i : orthogonal factor, alpha_i = <xp_{i-1}, p_{i-1}>

% normalized pi_1 = 1 before diving into iteration
beta = sqrt(polyint(pi_1, pi_1));
pi_1 = pi_1 / beta;
for i = 1:n
  alpha = polyint(conv(x, pi_1), pi_1);
  T(i, i) = alpha;
  pi = polymulti(x, pi_1);
  pi = polyminus(pi, alpha * pi_1);
  pi = polyminus(pi, beta*pi_2);
  beta = sqrt(polyint(pi, pi));
  pi = pi / beta;
  if (sum(abs(pi)) < TOL_ZERO | i ~= n)
    T(i, i+1) = beta;
    T(i+1, i) = beta;
    pi_2 = pi_1; pi_1 = pi;
  end
end


[S, Theta] = schur(T);
theta = (diag(Theta))';
s = S(1, :).^2;
s = s*2;
% Then, Gauss quadrature rule with n nodes becomes
%     \int_0^1 f(x) = \sum_{i=1}^n S(1, i) * f(theta_i),
% which should produce exact result for all polynomials
% with degree no more than 2*n-1.
%
% We test it here.

% ptest = zeros(1, 2*n-1);
% for itest = 1:2*n
% %  ptest = rand(2*n, 1);
%   ptest = zeros(1, 2*n);
%   ptest(2*n-itest+1) = 1;
%   intpexact = polyint(ptest, 1);
%   intpgauss = 0;
%   for i = 1:n
%     intpgauss = intpgauss + s(i) * polyval(ptest, theta(i));
%   end
%   % fprintf("x^%d, exact = %12.6f, Gauss = %12.6f.\n", ...
%   %         itest-1, intpexact, intpgauss);

%   if (abs(intpexact - intpgauss) > 1e-14)
%     fprintf("Wrong");
%     fprintf("x^%d, exact = %12.6f, Gauss = %12.6f.\n", ...
%             itest-1, intpexact, intpgauss);
%     pause(2);
%   end
% end

end

function out = polyint(p1, p2)
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

function poly = polyadd(p1, p2)
  p1 = reshape(p1, 1, []);
  p2 = reshape(p2, 1, []);
  degp1 = length(p1);
  degp2 = length(p2);
  if (degp1 < degp2)
    p1 = [zeros(1, degp2 - degp1), p1];
  else
    p2 = [zeros(1, degp1 - degp2), p2];
  end

  poly = p1 + p2;
end

function poly = polyminus(p1, p2)
  p1 = reshape(p1, 1, []);
  p2 = reshape(p2, 1, []);
  degp1 = length(p1);
  degp2 = length(p2);
  if (degp1 < degp2)
    p1 = [zeros(1, degp2 - degp1), p1];
  else
    p2 = [zeros(1, degp1 - degp2), p2];
  end

  poly = p1 - p2;
end

function poly = polymulti(p1, p2)
  poly = conv(p1, p2);
end
