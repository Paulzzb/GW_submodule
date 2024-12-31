function [results, time, relError, iter] = ...
      COmegaCstar_(Phi, Psi, evOcc, evUnocc, options)
% Calculate PhiPsi * Omega^{-1} * PhiPsi' by Cauchy integral, 
% where Omega_{ij} = evOcc_i - evUnocc_j, and PhiPsi = prod_states_gw(Phi, Psi);
% Inputs:
% Phi: wavefunctions, n * n1.  evOcc: n1.
% Psi: wavefunctions, n * n2.  evUnocc: n2.
% evOcc: all occupied states energies
% evUnocc: all unoccupied states energies
% options
% froErr: Tol for Cauchy integral 
% MaxIter: MaxIter for Cauchy integral. 
% Output:
% results: integral result.

% MaxIter = 12;
froErr = options.froErr;
MaxIter = options.MaxIter;
Error_array = [0, 0];
ii = sqrt(-1);
n1 = length(evOcc); nv = length(evUnocc);
n = size(Phi, 1);
M = evUnocc(end) - evOcc(end);
m = evUnocc(1) - evOcc(end);
if m <= 0
  error('No band gap in current system, Cauchy integral not applicable.')
end
k = (sqrt(M/m) - 1)/(sqrt(M/m)+1);
[K,Kstar] = ellipk(k);


results = zeros(n, n);
discretePoints = [-K, K];
newDiscretePoints = discretePoints;
Iker = @(s) 1 ./ sqrt((1+s.^2).*(1+k^2.*s.^2));
I = 1/2 * (integral(Iker,0,k^(-1)));

startTime = tic;

for iter = 1:MaxIter
  oldResults = results;
  results = results/2;
  
  % Numerical integral.
  for t = newDiscretePoints + I * ii
    [lambda, dlambda] = integrand(t, k, m, M);
    lambda = lambda + evOcc(end);
    OmegaOcc = diag(1./(lambda - evOcc));
    OmegaUnocc = diag(1./(lambda - evUnocc));
    OccMatrix = (Phi) * OmegaOcc * (Phi)';
    UnoccMatrix = conj(Psi) * OmegaUnocc * conj(Psi)';
%     OccMatrix = conj(Phi) * OmegaOcc * conj(Phi)';
%     UnoccMatrix = (Psi) * OmegaUnocc * (Psi)';
    matrix = (OccMatrix .* UnoccMatrix) .* (dlambda / pi /ii) ;
  
    % renew numerical integral, calculated by trapezoid formula.
    if iter == 1
      results = results + K*real(matrix);
    else
      results = results + newGap*real(matrix);
    end  
  end

  % Extrapolation.
  if iter >= 3
    oldResults_ = results_;
  end
  if iter >= 2
    results_ = (4*results - oldResults) ./ 3;
  end
  
  % Check exiting condition
  if iter >= 3
    Error = norm(oldResults_ - results_, 'fro') / norm(results_, 'fro');
    Error_array(iter) = Error;
    if Error <= froErr
      time = toc(startTime);
      results = results_;
      relError = Error;
      save('Cauchyintegralerror.mat', 'Error_array');
      return
    end
  end
  
  % Renew discreting points
  newGap = (discretePoints(2) - discretePoints(1)) / 2;
  newDiscretePoints = discretePoints(1:end-1) + newGap;
  oldDiscretePoints = discretePoints;
  discretePoints = zeros(1, 2^iter+1);
  discretePoints(1:2:end) = oldDiscretePoints;
  discretePoints(2:2:end) = newDiscretePoints;  
end

% Reach MaxIter, Print warning message.
save('Cauchyintegralerror.mat', 'Error_array');
warning(['Cauchy integral not converged after',num2str(MaxIter), 'iterations.']);
fprintf('Current relative error = %.5e, threshold = %.2e', Error, froErr);
relError = Error;
results = results_;
time = toc(startTime);
end % end of the function COmegaCstar



function [lambda, dlambda] = integrand(t, k, m, M)
% Calculate lambda and dlambda in Lu2017ISDFTPA
L = -log(k)/pi;
[SN,CN,DN] = ellipjc(t,L); 
lambda = sqrt(m*M).*((k^(-1)+SN)./(k^(-1)-SN));
dlambda = CN.*DN.*sqrt(m*M).*((2*k.^(-1))./(k^(-1)-SN).^2);
end % end of function integrand


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions below is copied directly from toolbox.


function [sn,cn,dn] = ellipjc(u,L,flag)
%ELLIPJC Jacobi elliptic functions for complex argument.
%   [SN,CN,DN] = ELLIPJC(U,L) returns the values of the Jacobi
%   elliptic functions evaluated at complex argument U and
%   parameter M=exp(-2*pi*L), 0 < L < Inf.  Recall that M = k^2,
%   where k is the elliptic modulus.
%   
%   U may be a matrix; L must be a scalar.  The entries of U are
%   expected to lie within the rectangle |Re U| < K, 0 < Im U <
%   Kp, where [K,Kp] = ELLIPK(L).
%   
%   Copyright (c) 1999 by Toby Driscoll. 
%   $Id: ellipjc.m 298 2009-09-15 14:36:37Z driscoll $

%   The built-in ELLIPJ can't handle compelx arguments, and
%   standard transformations to handle this would require ELLIPJ
%   called with parameter 1-M. When M < eps (or is even close),
%   this can't be done accurately.
%   
%   The algorithm is the descending Landen transformation,
%   described in L. Howell's PhD thesis from MIT. Additional
%   formulas from Gradshteyn & Ryzhik, 5th ed., and Abramowitz
%   & Stegun.

if nargin < 3
  % Absence of flag parameter indicates we must check for and transform u in
  % the upper half of the rectangle.
  [K,Kp] = ellipkkp(L);
  high = imag(u) > Kp/2;
  u(high) = i*Kp - u(high);
  m = exp(-2*pi*L);
else
  % Recursive call--L is actually m.
  high = zeros(size(u));
  m = L;
end

if m < 4*eps
  sinu = sin(u);
  cosu = cos(u);
  sn = sinu + m/4*(sinu.*cosu-u).*cosu;
  cn = cosu + m/4*(-sinu.*cosu+u).*sinu;
  dn = 1 + m/4*(cosu.^2-sinu.^2-1);
else
  if m > 1e-3
    kappa = (1-sqrt(1-m))/(1+sqrt(1-m));
  else
    kappa = polyval([132,42,14,5,2,1,0],m/4);
  end
  mu = kappa^2;
  v = u/(1+kappa);
  [sn1,cn1,dn1] = ellipjc(v,mu,1);
  denom = (1+kappa*sn1.^2);
  sn = (1+kappa)*sn1 ./ denom;
  cn = cn1.*dn1 ./ denom;
  dn = (1-kappa*sn1.^2) ./ denom;
end

if any(high(:))
  snh = sn(high);
  cnh = cn(high);
  dnh = dn(high);
  sn(high) = -1./(sqrt(m)*snh);
  cn(high) = i*dnh./(sqrt(m)*snh);
  dn(high) = i*cnh./snh;
end

end % end of function ellipjc



function [K,Kp] = ellipkkp(L)
%ELLIPKKP Complete elliptic integral of the first kind, with complement.
%   K = ELLIPKKP(L) returns the value of the complete elliptic
%   integral of the first kind, evaluated at M=exp(-2*pi*L), 0 < L
%   < Inf.
%   
%   [K,KP] = ELLIPKKP(L) also returns the result for complementary
%   parameter 1-M, which is useful when M < EPS.  Even when M <
%   1e-6, the built-in ELLIPKE can lose digits of accuracy for KP.
% 
%   Recall that the elliptic modulus k is related to the parameter
%   M by M = k^2.
% 
%   Copyright (c)1999 by Toby Driscoll. 
%   $Id: ellipkkp.m 298 2009-09-15 14:36:37Z driscoll $ 

%   ELLIPKKP uses the method of the arithmetic-geometric mean described
%   in 17.6 of M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%   Functions," Dover, 1965.  Same method as in ELLIPKE, only
%   interchanging 1 and 1-m to find KP.

% When m=exp(-2*pi*L) is extremely small, use O(m) approximations.
if L > 10
  K = pi/2;
  Kp = pi*L + log(4);
  return
end

m = exp(-2*pi*L);
a0 = 1;
b0 = sqrt(1-m);
s0 = m;
i1 = 0; mm = 1;
while mm > eps
  a1 = (a0+b0)/2;
  b1 = sqrt(a0.*b0);
  c1 = (a0-b0)/2;
  i1 = i1 + 1;
  w1 = 2^i1*c1.^2;
  mm = max(max(w1));
  s0 = s0 + w1;
  a0 = a1;
  b0 = b1;
end
K = pi./(2*a1);

im = find(m==1);
if ~isempty(im)
  K(im) = K(im)*inf;
end

if nargout > 1
  a0 = 1;
  b0 = sqrt(m);
  s0 = 1-m;
  i1 = 0; mm = 1;
  while mm > eps
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);
    c1 = (a0-b0)/2;
    i1 = i1 + 1;
    w1 = 2^i1*c1.^2;
    mm = max(max(w1));
    s0 = s0 + w1;
    a0 = a1;
    b0 = b1;
  end
  Kp = pi./(2*a1);
  im = find(m==0);
  if ~isempty(im)
    Kp(im) = Kp(im)*inf;
  end
end

end % end of ellipkkp
