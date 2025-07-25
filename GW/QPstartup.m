function QPstartup()
%   GW_STARTUP  Startup file for GW calculation
%   MAKE adds paths of the GW submodule to Matlab and choose a version to compute.


% CPATH = mfilename('fullpath');
% CPATH = fileparts(CPATH);
% CPATH = [CPATH, '/'];
restoredefaultpath;

% Set package root path
setappdata(0, 'PackageRoot', fileparts(mfilename('fullpath')));
CPATH = [getappdata(0, 'PackageRoot'), '/'];
disp(['GW module root path set to: ', CPATH]);

add_mpaths_only([CPATH 'common/']);
add_mpaths_only([CPATH 'driver_profile/']);
add_mpaths_only([CPATH 'input/']);
% add_mpaths_only([CPATH 'src_profile/']);
addpath(genpath([CPATH 'src_profile/']));
add_mpaths_only([CPATH 'GW_profile/']);
add_mpaths_only([CPATH 'test_profile/']);
add_mpaths_only([CPATH 'util_profile/']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Developer Hook] Insert your custom folders below
% add_mpaths_only([CPATH 'mymodule/']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Choose a version (default CPU)
% version = switchver('CPU');

end

function add_mpaths_only(pathroot)
% Add only directories that contain .m files (skip .mat/.data folders)

    % Get all subfolders
    all_dirs = strsplit(genpath(pathroot), pathsep);

    for i = 1:length(all_dirs)
        d = all_dirs{i};
        if isempty(d), continue; end

        % List .m files
        files = dir(fullfile(d, '*.m'));
        if ~isempty(files)
            addpath(d);
        end
    end
end
