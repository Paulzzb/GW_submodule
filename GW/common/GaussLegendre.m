function [weight, nodes] = GaussLegendre(n, tol)
%     [weight, nodes] = GaussLegendre(n, tol)
% Generate the Gauss quadrature rule on [-1, 1]
% Description:
% Gauss--Legendre quadrature rule generate a set of quadrature weights
% {c_k} and quadrature nodes {x_k}, and the integral can be approximated
% with 
% \int_{-1}^1 f(x) dx = \sum_{k=1}^n c_k f(x_k) + r_n(f),
% where the residual r_n(f) is bound by 
% |r_n(f)| <= f^(2n)(x)/(2*n!) \int_{-1}^1 \multi_{k=1}^n (x-x_k)^2 dx. 
% 
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
%     weight: weights
%     nodes: nodes

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

clist = zeros(n, 1);
alist = zeros(n, 1);
for i = 1:n
  alist(i) = (2*i-1) / i;
  clist(i) = (i-1) / i;
end

J = zeros(n);
for i = 1:n-1
  beta = sqrt(clist(i+1) / alist(i) / alist(i+1));
  J(i, i+1) = beta;
  J(i+1, i) = beta;
end

[qlist, tlist] = eig(J);
tlist = diag(tlist);
nodes = tlist';
weight = qlist(1, :).^2 * 2;

end % EOF



