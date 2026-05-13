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
add_mpaths_only([CPATH 'tmp_profile/']);
% add_mpaths_only([CPATH 'service/']);
addpath(genpath([CPATH 'service/']))
addpath(genpath([CPATH 'packages/']))
% add_mpaths_only([CPATH 'database_profile/']);
addpath(genpath([CPATH 'database_profile/']))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Developer Hook] Insert your custom folders below
% add_mpaths_only([CPATH 'mymodule/']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optional ISDF MEX: compile if missing/outdated, then report availability.
qp_setup_isdf_mex(CPATH);

% Choose a version (default CPU)
% version = switchver('CPU');

end

function qp_setup_isdf_mex(CPATH)
%QP_SETUP_ISDF_MEX Build and verify ISDF-related MEX (prod + Schur rank-1).
%
% Set environment variable QP_SKIP_ISDF_MEX to any non-empty value to skip
% compilation attempts (verification is still printed).

  if ~isempty(getenv('QP_SKIP_ISDF_MEX'))
    fprintf('[QPstartup] ISDF MEX: build skipped (QP_SKIP_ISDF_MEX is set).\n');
    qp_report_isdf_mex();
    return;
  end

  ext = ['.', mexext];
  isdfDir = fullfile(CPATH, 'service', '+isdf');
  adaptiveDir = fullfile(isdfDir, '+adaptive');

  targets = {
    fullfile(isdfDir, 'prod_mex.c'), fullfile(isdfDir, ['prod_mex', ext]), @isdf.build_prod_mex;
    fullfile(adaptiveDir, 'isdf_schur_rank1_mex.c'), fullfile(adaptiveDir, ['isdf_schur_rank1_mex', ext]), @isdf.adaptive.build_schur_rank1_mex;
    fullfile(adaptiveDir, 'isdf_schur_rank1_prod_mex.c'), fullfile(adaptiveDir, ['isdf_schur_rank1_prod_mex', ext]), @isdf.adaptive.build_schur_rank1_prod_mex
    };

  for k = 1:size(targets, 1)
    cfile = targets{k, 1};
    mexfile = targets{k, 2};
    buildFcn = targets{k, 3};
    if ~exist(cfile, 'file')
      warning('QPstartup:IsdfMexSrc', 'Missing C source, skip: %s', cfile);
      continue;
    end
    if qp_mex_needs_rebuild(cfile, mexfile)
      try
        fprintf('[QPstartup] Building MEX: %s\n', mexfile);
        buildFcn();
      catch ME
        warning('QPstartup:IsdfMexBuild', 'MEX build failed for %s: %s', mexfile, ME.message);
      end
    end
  end

  qp_report_isdf_mex();
end

function tf = qp_mex_needs_rebuild(cfile, mexfile)
  if ~exist(cfile, 'file')
    tf = false;
    return;
  end
  if ~exist(mexfile, 'file')
    tf = true;
    return;
  end
  dc = dir(cfile);
  dm = dir(mexfile);
  tf = datenum(dc(1).date) > datenum(dm(1).date);
end

function qp_report_isdf_mex()
  rows = {
    'isdf.prod_mex';
    'isdf.adaptive.isdf_schur_rank1_mex';
    'isdf.adaptive.isdf_schur_rank1_prod_mex'
    };
  fprintf('[QPstartup] ISDF MEX availability (which):\n');
  for i = 1:numel(rows)
    name = rows{i};
    w = which(name);
    if isempty(w)
      fprintf('  %-42s %s\n', name, '(missing)');
    else
      fprintf('  %-42s %s\n', name, w);
    end
  end
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
