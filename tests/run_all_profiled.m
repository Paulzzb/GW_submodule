% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/10 ZZ

function run_all_profiled(out_dir)
%RUN_ALL_PROFILED  Run example/run_all under MATLAB profiler and save HTML.
%
%   run_all_profiled
%   run_all_profiled(out_dir)
%
%   Default out_dir: example/profile_run_all/
%   Writes profiler HTML via profsave (open out_dir/file0.html).

  example_dir = fileparts(mfilename('fullpath'));
  if isempty(example_dir)
    example_dir = pwd;
  end
  addpath(example_dir);

  if nargin < 1 || isempty(out_dir)
    out_dir = fullfile(example_dir, 'profile_run_all');
  end
  if ~exist(out_dir, 'dir')
    mkdir(out_dir);
  end

  fprintf('=== run_all_profiled ===\n');
  fprintf('profile out = %s\n', out_dir);

  profile('on');
  try
    run_all;
  catch ME
    info = profile('info');
    profile('off');
    try
      profsave(info, out_dir);
    catch
    end
    rethrow(ME);
  end

  info = profile('info');
  profile('off');
  profsave(info, out_dir);
  fprintf('Profile saved: %s (open file0.html)\n', out_dir);
end
