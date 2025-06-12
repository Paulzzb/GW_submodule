function GWlog(msg, level)
% GWlog Unified logging function for GW toolbox
%   GWlog(MSG, LEVEL)
%   LEVEL = 0: always shown
%   LEVEL = 1: normal progress messages
%   LEVEL = 2: verbose/debug messages

    persistent VERBOSE
    if isempty(VERBOSE)
        VERBOSE = 1; % Default verbosity level
    end

    if level <= VERBOSE
        fprintf('[GW] %s\n', msg);
    end
end
