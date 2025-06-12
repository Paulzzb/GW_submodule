function gw_startup()
%   GW_STARTUP  Startup file for GW calculation
%   MAKE adds paths of the GW submodule to Matlab and choose a version to compute.


% CPATH = mfilename('fullpath');
% CPATH = fileparts(CPATH);
% CPATH = [CPATH, '/'];
restoredefaultpath;

setappdata(0, 'PackageRoot', fileparts(mfilename('fullpath')));

disp(['GW module root path set to: ', getappdata(0, 'PackageRoot')]);

CPATH = [getappdata(0, 'PackageRoot'), '/'];

addpath(genpath([CPATH 'common/']));

addpath(genpath([CPATH 'driver_profile/']));

addpath(genpath([CPATH 'input/']));

addpath(genpath([CPATH 'src_profile/']));

addpath(genpath([CPATH 'GW_profile/']));

addpath(genpath([CPATH 'test_profile/']));

addpath(genpath([CPATH 'util_profile/']));


% Choose a version (default CPU)
% version = switchver('CPU');

end
