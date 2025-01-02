function KSSOLV_startup()
% KSSOLV_STARTUP  Startup file for KSSOLV
%   MAKE adds paths of the KSSOLV to Matlab and choose a version to compute.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

file_path = mfilename('fullpath');
tmp = strfind(file_path,'KSSOLV_startup');
file_path = file_path(1:(tmp(end)-1));

% Folder for all utility functions
addpath([file_path 'util']);

% Foulder for all source files recursively
addpath(genpath([file_path 'src']));

% Foulder for all external files recursively
addpath(genpath([file_path 'external']));

% Foulder for all source files recursively
addpath(genpath([file_path 'test']));

addpath(genpath([file_path 'structure']));

% Choose a version (default CPU)
version = switchver('CPU');

end
