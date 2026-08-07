% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/07/27 ZZ

function QPstartup()
%QPSTARTUP Startup for the GW module.
%   Add package paths and ensure required ISDF MEX kernels are available.

% Reset MATLAB path to a clean state before adding package folders.
restoredefaultpath;

% Record package root for later use by drivers and services.
setappdata(0, 'PackageRoot', fileparts(mfilename('fullpath')));
CPATH = [getappdata(0, 'PackageRoot'), '/'];
disp(['GW module root path set to: ', CPATH]);

% Add core folders to the MATLAB path.
addpath(CPATH)
add_mpaths_only([CPATH 'common/']);
add_mpaths_only([CPATH 'driver/']);
add_mpaths_only([CPATH 'input/']);
addpath(genpath([CPATH 'src/']));
add_mpaths_only([CPATH 'util/']);
add_mpaths_only([CPATH 'example/']);
add_mpaths_only([CPATH 'test_profile/']);
addpath(genpath([CPATH 'service/']))
addpath(genpath([CPATH 'packages/']))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Developer Hook] Insert your custom folders below
% add_mpaths_only([CPATH 'mymodule/']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build missing ISDF MEX kernels for the current platform after paths are ready.
% which() only sees mex for this OS; if absent, compile via mex() into the package dir.
mex_specs = {
  'isdf.prod_mex', @() isdf.build_prod_mex
  'isdf.adaptive_double.isdf_schur_rank1_mex', @() isdf.adaptive_double.build_schur_rank1_mex
  'isdf.adaptive_double.isdf_schur_rank1_prod_mex', @() isdf.adaptive_double.build_schur_rank1_prod_mex
  'isdf.adaptive_single.isdf_schur_rank1_mex', @() isdf.adaptive_single.build_schur_rank1_mex
  'isdf.adaptive_single.isdf_schur_rank1_prod_mex', @() isdf.adaptive_single.build_schur_rank1_prod_mex
  };
for i = 1:size(mex_specs, 1)
  ensure_mex_kernel(mex_specs{i, 1}, mex_specs{i, 2});
end

% Report MEX availability for quick diagnostics.
fprintf('MEX availability (platform=%s):\n', computer('arch'));
for i = 1:size(mex_specs, 1)
  sym = mex_specs{i, 1};
  fprintf('  %-52s %d\n', sym, ~isempty(which(sym)));
end

end

function add_mpaths_only(pathroot)
%ADDPATHS_ONLY Add only subfolders that contain .m files.

    all_dirs = strsplit(genpath(pathroot), pathsep);

    for i = 1:length(all_dirs)
        d = all_dirs{i};
        if isempty(d), continue; end

        files = dir(fullfile(d, '*.m'));
        if ~isempty(files)
            addpath(d);
        end
    end
end

function ensure_mex_kernel(symbol_name, build_fn)
%ENSURE_MEX_KERNEL Build a MEX kernel if it is not on the MATLAB path.
%   Detects the current platform binary via which(); compiles with mex() if missing.

  if ~isempty(which(symbol_name))
    return;
  end

  fprintf('Building missing MEX: %s\n', symbol_name);
  try
    build_fn();
    rehash path;
  catch ME
    warning('QPstartup:MexBuildFailed', ...
      'Failed to build %s: %s', symbol_name, ME.message);
    return;
  end

  if isempty(which(symbol_name))
    warning('QPstartup:MexStillMissing', ...
      'Built %s but which() still empty (check mex setup / path).', symbol_name);
  end
end

