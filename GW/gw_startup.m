function gw_startup()
%   GW_STARTUP  Startup file for GW calculation
%   MAKE adds paths of the GW submodule to Matlab and choose a version to compute.


% CPATH = mfilename('fullpath');
% CPATH = fileparts(CPATH);
% CPATH = [CPATH, '/'];

setappdata(0, 'PackageRoot', fileparts(mfilename('fullpath')));

disp(['GW module root path set to: ', getappdata(0, 'PackageRoot')]);

CPATH = [getappdata(0, 'PackageRoot'), '/'];



% Folder for all utility functions
addpath(genpath([CPATH 'util/']));

% Foulder for all source files recursively
addpath(genpath([CPATH 'src/']));

% Foulder for all external files recursively
addpath(genpath([CPATH 'test/']));

% Foulder for all source files recursively
addpath(genpath([CPATH 'example/']));

% Foulder for all source files recursively
addpath(genpath([CPATH 'data/']));
% addpath(genpath([CPATH 'structure']));

% Choose a version (default CPU)
% version = switchver('CPU');

end
