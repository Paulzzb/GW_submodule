function [x, f, g, frec, grec, alpharec] = nlcg(pars, options)
% NLCG Nonlinear Conjugate Gradient minimization.
%   Calls:  [x, f, g] = nlcg(pars) 
%    and:   [x, f, g, frec, alpharec] = nlcg(pars, options)
% 
%   Input parameters
%    pars is a required struct, with two required fields
%      pars.nvar: the number of variables
%      pars.fgname: string giving the name of function (in single quotes) 
%         that returns the function and its gradient at a given input x, 
%         with call   [f,g] = fgtest(x,pars)  if pars.fgname is 'fgtest'.
%         Any data required to compute the function and gradient may be
%         encoded in other fields of pars.
%    options is an optional struct, with no required fields
%      options.x0: each column is a starting vector of variables
%          (default: empty)
%      options.nstart: number of starting vectors, generated randomly
%          if needed to augment those specified in options.x0
%          (default: 10 if options.x0 is not specified)
%      options.maxit: max number of iterations 
%          (default 1000) (applies to each starting vector)
%      options.normtol: termination tolerance on gradient norm 
%          (default: 1e-6) (applies to each starting vector)
%      options.fvalquit: quit if f drops below this value 
%          (default: -inf) (applies to each starting vector)
%      options.xnormquit: quit if norm(x) exceeds this value
%          (default: inf)
%      options.cpumax: quit if cpu time in secs exceeds this 
%          (default: inf) (applies to total running time)
%      options.version:
%           'P' for Polak-Ribiere-Polyak (not recommended: fails on hard problems)
%           'F' for Fletcher-Reeves (not recommended: often stagnates)
%           'C' for Polak-Ribiere-Polyak Constrained by Fletcher-Reeves 
%              (recommended, combines advantages of 'P' and 'F'; default)
%           'S' for Hestenes-Stiefel (not recommended)
%           'Y' for Dai-Yuan (allows weak Wolfe line search, see below)
%           'Z' for Hager-Zhang
%           '-' for Steepest Descent (for comparison)
%      options.strongwolfe: 1 for strong Wolfe line search (default)
%                           0 for weak Wolfe line search (not recommended)
%      options.wolfe1: first Wolfe line search parameter 
%          (ensuring sufficient decrease in function value, default: 1e-4)
%      options.wolfe2: second Wolfe line search parameter 
%          (strong Wolfe condition: ensuring decrease in absolute value of
%           projected gradient; must be < 1/2 to guarantee convergence;
%           default: 0.49)
%          ("weak" Wolfe line search is not recommended for use with
%           most versions; the convergence theory requires a "strong" Wolfe
%           search, but version 'Y' (Dai-Yuan) is an exception)
%      options.prtlevel: one of 0 (no printing), 1 (minimal), 2 (verbose)
%          (default: 1)
%
%   Output parameters 
%    x: each column is an approximate minimizer, one for each starting vector
%    f: each entry is the function value for the corresponding column of x
%    g: each column is the gradient at the corresponding column of x
%    frec: f values for all iterates (not including those in the line search)
%     (cell array, frec{j} is for jth starting vector)
%    alpharec: record of steps taken in the line search 
%     (also cell array)
%
%   NLCG (Nonlinear Conjugate gradient)
%   is intended for minimizing smooth, not necessarily convex, functions.
%   The Fletcher-Reeves version (version='F') is globally convergent in
%   theory but often stagnates in practice.  The Polak-Ribiere version
%   (version='P') works better in practice but its search direction
%   may not even be a descent direction and it may not converge.  
%   The 'C' version combines the best of both.  It is Polak-Ribiere 
%   constrained by Fletcher-Reeves, typically behaving as well as or better
%   than PR in practice but with the same global convergence guarantee as FR.
%   The Hestenes-Stiefel version (version='S') also often fails.
%   The Dai-Yuan (version='Y') is the first to allow use of a weak Wolfe
%   line search.  The Hager-Zhang (version='Z') is the newest and is also
%   promising.
%   Reference: second edition of Nocedal and Wright, Chapter 5, plus papers
%   by Dai-Yuan and Hager-Zhang.
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "nlcg".
%   NLCG Version 1.0, 2010, see GPL license info below.

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

% parameter defaults
if nargin == 0
   error('pars is a required input parameter')
end
if nargin == 1
   options = [];
end
options = setdefaults(pars,options);  % set default options
options = setx0(pars, options); % augment options.x0 randomly
cpufinish = cputime + options.cpumax;
fvalquit = options.fvalquit;
xnormquit = options.xnormquit;
prtlevel = options.prtlevel;
% set CG options, including version and line search options
options = setdefaultsnlcg(pars, options); 
x0 = options.x0;
nstart = size(x0,2);
for run = 1:nstart
    if prtlevel > 0 && nstart > 1
        fprintf('nlcg: starting point %d\n', run)
    end
    options.cpumax = cpufinish - cputime; % time left
    [x(:,run), f(run), g(:,run), frec{run}, grec{run}, alpharec{run}] = nlcg1run(x0(:,run), pars, options);
    if cputime > cpufinish || f < fvalquit || norm(x) > xnormquit
        break
    end
end