function options = setdefaultsnlcg(pars, options)
%  set defaults for NLCG(in addition to those already set by setdefaults)
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "nlcg".
%   NLCG Version 1.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  NLCG 2.0 Copyright (C) 2010  Michael Overton
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
% line search options
if isfield(options, 'strongwolfe')
    if options.strongwolfe ~= 0 & options.strongwolfe ~= 1
        error('setdefaultsnlcg: input "options.strongwolfe" must be 0 or 1')
    end
else
    % strong Wolfe is very complicated but is normally needed for NLCG
    options.strongwolfe = 1;  
end
if isfield(options, 'wolfe1') % conventionally anything in (0,1), but 0 is OK while close to 1 is not
    if ~isreal(options.wolfe1) | options.wolfe1 < 0 | options.wolfe1 > 0.5
        error('setdefaultsnlcg: input "options.wolfe1" must be between 0 and 0.5')
    end
else
    options.wolfe1 = 0; % conventionally this should be positive, but zero is fine in practice
end
if isfield(options, 'wolfe2') % conventionally should be > wolfe1, but both 0 OK for e.g. Shor
    if ~isreal(options.wolfe2) | options.wolfe2 < options.wolfe1  | options.wolfe2 >= 1
        error('setdefaultsnlcg: input "options.wolfe2" must be between max(0,options.wolfe1) and 1')
    end
else % this can be set to zero to simulate an exact line search
    options.wolfe2 = 0.5;  % should be < 1/2 for convergence theory for FR, PRFR
end
if options.strongwolfe == 0
    if ~exist('linesch_ww')
        error('"linesch_ww" is not in path: it can be obtained from the HANSO distribution')
    end
else
    if ~exist('linesch_sw')
        error('"linesch_sw" is not in path: it is required for strong Wolfe line search')
    end
end
if isfield(options, 'quitLSfail')
    if options.quitLSfail ~= 0 & options.quitLSfail ~= 1
        error('setdefaultsnlcg: input "options.quitLSfail" must be 0 or 1')
    end
else
    if options.strongwolfe == 1 & options.wolfe2 == 0
        % simulated exact line search, so don't quit if it fails
        options.quitLSfail = 0;
    else  % quit if line search fails
        options.quitLSfail = 1;
    end
end
% which version of CG
if isfield(options, 'version')
    version = options.version;
else 
    version = 'C'; % default is Polak-Ribiere constrained by Fletcher-Reeves
end
version = upper(version);
if version == 'F' | version == 'P' | version == 'C' | ...
  version == 'S'| version == 'Y' | version == 'Z' | version == '-'
    % nothing to do
else
    error('nlcg: options.version must be F, P, C, S, Y, Z or - (use single quotes)')
end
options.version = version;
prtlevel = options.prtlevel;
if prtlevel > 0
    if version ~= '-'
        fprintf('Conjugate Gradient, version %s', version)
    else
        fprintf('Steepest Descent')
    end
end
if prtlevel > 0
    if options.strongwolfe
        if options.wolfe2 > 0
            fprintf(' with strong Wolfe line search\n')
        else
            fprintf(' simulating exact line search\n')
        end
    else
        fprintf(' with weak Wolfe line search (NOT generally recommended for CG) \n')
    end
end