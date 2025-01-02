% n = 5;
tol = 1e-8;

% Generate the Gauss quadrature rule on [0, 1]

% % A(i, j) = \int_0^1 x^(i-1)x^(j-1) dx = 1 / (i+j-1)
% A = zeros(n, n);
% for i = 1:n
%   for j = 1:n
%     A(i, j) = 1 / (i+j-1);
%   end
% end
%
% % Do cholesky decomposition to A
% R = chol(A, 'upper');
% S = inv(R);

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
  if abs(imag(beta)) > 1e-10
    ;
  end
  pi = pi / beta;
  if (sum(abs(pi)) < tol | i ~= n)
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

ptest = zeros(1, 2*n-1);
for itest = 1:2*n
%  ptest = rand(2*n, 1);
  ptest = zeros(1, 2*n);
  ptest(2*n-itest+1) = 1;
  intpexact = polyint(ptest, 1);
  intpgauss = 0;
  for i = 1:n
    intpgauss = intpgauss + s(i) * polyval(ptest, theta(i));
  end
  % fprintf("x^%d, exact = %12.6f, Gauss = %12.6f.\n", ...
  %         itest-1, intpexact, intpgauss);

  if (abs(intpexact - intpgauss) > 1e-14)
    % fprintf("Wrong");
    % fprintf("x^%d, exact = %12.6f, Gauss = %12.6f.\n", ...
            % itest-1, intpexact, intpgauss);
    % pause(2);
    ;
  end
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