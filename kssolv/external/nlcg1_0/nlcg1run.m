function [x, f, g, frec, grec, alpharec] = nlcg1run(x0, pars, options)
% make a single run of nonlinear conjugate gradient from one starting point
% intended to be called by nlcg only
% Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
% with a subject header containing the string "nlcg".
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

fgname = pars.fgname;
normtol = options.normtol;
fvalquit = options.fvalquit;
cpufinish = cputime + options.cpumax;
maxit = options.maxit;
prtlevel = options.prtlevel;
version = options.version;
strongwolfe = options.strongwolfe;
wolfe1 = options.wolfe1;
wolfe2 = options.wolfe2; 
quitLSfail = options.quitLSfail;
frec = nan; % so defined in case of immediate return
grec = nan;
alpharec = nan;
x = x0;
[f,g] = feval(fgname, x, pars);
if size(g,2) > size(g,1) % error return is appropriate here
    error('gradient must be returned as a column vector, not a row vector')
end
gnorm = norm(g);


if f == inf % better not to generate an error return
    if prtlevel > 0
        fprintf('nlcg: f is infinite at initial iterate\n')
    end
    return
elseif isnan(f)
    if prtlevel > 0
        fprintf('nlcg: f is nan at initial iterate\n')
    end
    return
elseif gnorm < normtol
    if prtlevel > 0
       fprintf('nlcg: tolerance on gradient satisfied at initial iterate\n')
    end
    return
elseif f < fvalquit
    if prtlevel > 0
        fprintf('nlcg: below target objective at initial iterate\n')
    end
    return
end 
p = -g;  % start with steepest descent

grec = [];                      % Added by Subhajit; 06/2017; LBNL

for iter = 1:maxit
    gtp = g'*p;
    if  gtp >= 0 || isnan(gtp)
        if prtlevel > 0
            fprintf('nlcg: not descent direction, quit at iteration %d, f = %g, gnorm = %5.1e\n',...
                iter, f, gnorm)
        end
        return
    end
    gprev = g;
    if strongwolfe % strong Wolfe line search is usually essential (default)
        [alpha, x, f, g, fail] = ...
            linesch_sw(x, f, g, p, pars, wolfe1, wolfe2, fvalquit, prtlevel);
    else  % still, leave weak Wolfe as an option for comparison
        [alpha, x, f, g, fail] = ...
            linesch_ww(x, f, g, p, pars, wolfe1, wolfe2, fvalquit, prtlevel);
    end
    gnorm = norm(g);
    frec(iter) = f;
    grec = [grec g];           % Added by Subhajit; 06/2017; LBNL
    % grec(:,iter) = g;        % Added by Subhajit; 06/2017; LBNL
    alpharec(iter) = alpha;
    
  
    fprintf('Forces for all atoms=');
    force = reshape(-g, [numel(g)/3, 3]);
    disp(force);

    if prtlevel > 1
        fprintf('nlcg: iter %d: step = %5.1e, f = %g, gnorm = %5.1e\n', iter, alpha, f, gnorm)
    end
    if f < fvalquit
        if prtlevel > 0
            fprintf('nlcg: reached target objective, quit at iteration %d \n', iter)
        end
        return
    end
    if fail == 1 % Wolfe conditions not both satisfied, quit
        if ~quitLSfail 
            if prtlevel > 1
                fprintf('nlcg: continue although line search failed\n')
            end
        else % quit since line search failed
            if prtlevel > 0
                fprintf('nlcg: line search conds not satisfied, quit at iteration %d, f = %g, gnorm = %5.1e\n',...
                  iter, f, gnorm)
            end
            return
        end
    elseif fail == -1 % function apparently unbounded below
        if prtlevel > 0
            fprintf('nlcg: function may be unbounded below, quit at iteration %d, f = %g\n', iter, f)
        end
        return
    end
    if gnorm <= normtol
        if prtlevel > 0
            fprintf('nlcg: gradient norm below tolerance, quit at iteration %d, f = %g\n', iter, f)
        end
        return
    end
    if cputime > cpufinish
        if prtlevel > 0
            fprintf('nlcg: cpu time limit exceeded, quit at iteration %d, f = %g\n', iter, f)
        end
        return
    end
    y = g - gprev;
    if version == 'P'
        nmgprevsq = gprev'*gprev;
        beta = (g'*y)/nmgprevsq;  % Polak-Ribiere-Polyak
    elseif version == 'F'
        nmgprevsq = gprev'*gprev;
        beta = (g'*g)/nmgprevsq;  % Fletcher-Reeves  
    elseif version == 'C' % combined PR-FR (suggested by Gilbert-Nocedal)
    % ensures beta <= |beta_fr|, allowing proof of % global convergence,
    % but avoids inefficiency of FR which happens when beta_fr gets stuck near 1
        nmgprevsq = gprev'*gprev;
        beta_pr = (g'*y)/nmgprevsq;  % Polak-Ribiere-Polyak
        beta_fr = (g'*g)/nmgprevsq;  % Fletcher-Reeves  
        if beta_pr < -beta_fr  
            if prtlevel > 1
                fprintf('*** truncating beta_pr = %g to -beta_fr = %g\n', beta_pr, -beta_fr)
            end
            beta = -beta_fr;  
        elseif beta_pr > beta_fr
            if prtlevel > 1
                fprintf('*** truncating beta_pr = %g to +beta_fr = %g\n', beta_pr, beta_fr)
            end
            beta = beta_fr;
        else
            beta = beta_pr;
        end
    elseif version == 'S' % Hestenes-Stiefel
        beta = (g'*y)/(p'*y);
    elseif version == 'Y' % Dai-Yuan
        beta = (g'*g)/(p'*y);
    elseif version == 'Z' % Hager-Zhang
        pty = p'*y; % p is called d in their paper
        theta = 2*(y'*y)/pty;
        beta_hz = ((y-theta*p)'*g)/pty;
        eta = -1/(norm(p)*min(.01,norm(gprev)));
        beta = max(beta_hz, eta);
    elseif version == '-' % Steepest Descent
        beta = 0;
    else
        error('nlcg: no such version') % already checked in nlcg
    end
    p = beta*p - g;

  
  end % for loop
if prtlevel > 0
    fprintf('nlcg: %d iterations reached, f = %g, gnorm = %5.1e\n', maxit, f, gnorm)
end
