function KSSOLV_startup()
% KSSOLV_STARTUP  Startup file for KSSOLV
%   MAKE adds paths of the KSSOLV to Matlab.

file_path = mfilename('fullpath');
tmp = strfind(file_path,'KSSOLV_startup');
file_path = file_path(1:(tmp(end)-1));

% Folder for all utility functions
%addpath([file_path 'kssolvsrc/Util/Tools']);
%addpath([file_path 'kssolvsrc/Util/Constants']);
%addpath([file_path 'kssolvsrc/Util/SpecialFunctions']);

% Folder for all pseuopotential files
addpath([file_path 'ppdata/default']);

% Folder for all source files recursively
addpath(genpath([file_path 'kssolvsrc']));

% Folder for all utility functions
%addpath([file_path 'kssolvsrc/Pseudopotential']);

% Folder for all external files recursively
addpath(genpath([file_path 'external']));
addpath(genpath([file_path 'examples']));

addpath(genpath([file_path 'gw']));

end
