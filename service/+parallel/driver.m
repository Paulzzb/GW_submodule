% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06

function driver(~, config)
  obj = parallel.base.parallel_m();
  obj.toolbox_available = license('test', 'Distrib_Computing_Toolbox');

  requested = true;
  workers = int32(0);

  if isstruct(config) && isfield(config, 'PARALLEL')
    pconf = config.PARALLEL;
    if isstruct(pconf)
      if isfield(pconf, 'enabled')
        requested = logical(pconf.enabled);
      end
      if isfield(pconf, 'workers')
        validateattributes(pconf.workers, {'numeric'}, ...
          {'scalar', 'integer', 'nonnegative'}, mfilename, 'config.PARALLEL.workers');
        workers = int32(pconf.workers);
      end
    end
  end

  obj.requested = requested;
  obj.workers = workers;
  obj.enabled = obj.toolbox_available && requested;
  if obj.enabled
    try
      pool = gcp('nocreate');
      target_workers = double(workers);
      if target_workers > 0
        if isempty(pool) || pool.NumWorkers ~= target_workers
          if ~isempty(pool)
            delete(pool);
          end
          c = parcluster('local');
          if c.NumWorkers < target_workers
            c.NumWorkers = target_workers;
          end
          parpool(c, target_workers);
        else
          fprintf('[parallel.driver] reuse existing pool (%d workers)\n', pool.NumWorkers);
        end
      else
        if isempty(pool)
          parpool('local');
        end
      end
      obj.mode = "parfor";
    catch ME
      warning('parallel:driver:pool', ...
        'Failed to initialize parallel pool; fall back to serial mode. Reason: %s', ME.message);
      obj.enabled = false;
      obj.mode = "serial";
    end
  else
    obj.mode = "serial";
  end

  parallel.save2mod(obj);
end
