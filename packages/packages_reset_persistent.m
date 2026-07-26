% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/24

function report = packages_reset_persistent()
  % Reset persistent caches across package functions.

  cleanup_list = {
    @() SCATTER_Bamp('reset'), 'SCATTER_Bamp(''reset'')';
  };

  report.ok = true;
  report.cleared = {};
  report.errors = {};

  for i = 1:size(cleanup_list, 1)
    fn = cleanup_list{i, 1};
    name = cleanup_list{i, 2};

    try
      fn();
      report.cleared{end + 1} = name; %#ok<AGROW>
    catch ME
      report.ok = false;
      report.errors{end + 1} = sprintf('%s: %s', name, ME.message); %#ok<AGROW>
    end
  end
end
