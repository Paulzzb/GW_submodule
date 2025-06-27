function [alpha, x, f, grad, fail, nsteps] = ...
           linesch_sw(x0, f0, grad0, d, pars, c1, c2, fvalquit, prtlevel);
%  based on Nocedal-Wright line search, see their book
%  strong Wolfe line search with cubic interpolation
%  not recommended for BFGS, especially if function may be nonsmooth: 
%    use the far simpler and preferable linesch_ww instead
%  recommended for use with CG, where strong Wolfe condition needed for
%  convergence analysis
% call:  [alpha, x, f, grad, fail, nsteps] = ...
%          linesch_sw(x0, f0, grad0, d, pars, c1, c2, prtlevel);
%  written by M. Overton (overton@cs.nyu.edu)
%  input
%   x0:      intial point
%   f0:      function value at x0
%   g0:      gradient at x0
%   d:       search direction (called p in Nocedal and Wright's book) 
%   pars:    a structure that specifies the function name as well
%            anything else that the user needs to access in programming the
%            function and gradient values
%   fgname:  name of function that returns function and gradient
%            it expects as input only x and pars, a parameter structure 
%            it is invoked by: [f,g] = feval(fgname, x, pars)
%   c1: Wolfe parameter for the condition 
%          f(x0 + alpha d) < = f0 + c1*alpha*grad0'*d   (default 0)
%        in theory, should be positive, but 0 works fine in practice
%   c2: Wolfe parameter for the (strong) condition 
%          abs((grad f)(x0 + alpha d)'*d) < = -c2*grad0'*d
%        these should satisfy 0 <= c1 <= c2 <= 1     (default 0.5)
%        usually, c1 should be small but c2 should not, but set c2 small
%         if an accurate line search is desired
%        for BFGS, only need c2 < 1 for update to be positive definite 
%        for nonlinear conjugate gradient, need c2 < 0.5 for the 
%         Fletcher-Reeves convergence analysis to apply
%   fvalquit: supposed to quit immediately if f drops below this value, 
%        regardless of the Wolfe conditions, to be consistent with
%        linesch_ww, but this is not implemented yet, so an error
%        is generated if fvalquit is not equal to -inf
%   prtlevel: 2 for info when line search fails (default), 0 for no printing
%      
%  output:
%   alpha: steplength satisfying Wolfe conditions
%   x:     x0 + alpha*d
%   f:     f(x0 + alpha d)
%   grad:  (grad f)(x0 + alpha d)
%   fail:  0 for normal termination (Wolfe conditions satisfied)
%          1 if too many steps were taken in lszoom (Wolfe conds not satisfied)
%          -1 if minimizer was not bracketed (too many expansion steps, 
%                                              function may be unbounded below)
%   nsteps: number of steps taken in lszoom
% NLCG Version 1.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  NLCG 1.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6  % check if the optional Wolfe parameters were passed
    c1 = 0;  % should be positive in theory, not necessary in practice
end
if nargin < 7
    c2 = 0.5;  % should be < 1/2 in theory
end
if nargin < 8
    fvalquit = -inf;
end
if nargin < 9
    prtlevel = 2;
end
if fvalquit ~= -inf
    warning('linesch_sw: fvalquit option has not been implemented')
end
if c1 < 0 | c1 > c2 | c2 > 1
   warning('linesch_sw: Wolfe parameters are supposed to satisfy 0 <= c1 <= c2 <= 1')
end
g0 = grad0'*d;
if g0 >= 0
   error('linesch_sw: g0 must be negative, not a descent direction')
end
old = 0; fold = f0; gold = g0;
new = 1;
nexpand = max([50  -round(log2(norm(d)))]);
fgname = pars.fgname;
for k = 1:nexpand
   xnew = x0 + new*d;
   [fnew,gradnew] = feval(fgname, xnew, pars);
   gnew = gradnew'*d;
   if fnew > f0 + c1*new*g0 | ((fnew >= fold) & k > 1)   % gone too far
      [alpha, x, f, grad, fail, nsteps] = lszoom(old, new, ...
          fold, fnew, gold, gnew, f0, g0, x0, d, pars, c1, c2, prtlevel);
      return
   end
   if abs(gnew) <= -c2*g0   % quit, conditions are satisfied first time
      alpha = new; x = xnew; f = fnew; grad = gradnew; fail = 0; nsteps = k;
      return
   end
   if gnew >= 0     % gone past minimizer, but new point is better than old
      [alpha, x, f, grad, fail, nsteps] = lszoom(new, old, ... %note reversal
          fnew, fold, gnew, gold, f0, g0, x0, d, pars, c1, c2, prtlevel);
      return
   end
   old = new;       % minimizer not bracketed yet
   fold = fnew;
   gold = gnew;
   new = 2*new;     
end % loop
if prtlevel > 1
    fprintf('linesch_sw: minimizer was not bracketed in %d expansion steps, ', nexpand)
    fprintf('function may be unbounded below\n')
end
alpha = new; 
x = xnew; 
f = fnew; 
grad = gradnew; 
fail = -1;
nsteps = 0;  % lszoom was never called