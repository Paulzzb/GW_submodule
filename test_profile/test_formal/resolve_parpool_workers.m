function n = resolve_parpool_workers()
%RESOLVE_PARPOOL_WORKERS  Worker count for parpool on this host / Slurm allocation.
%
%   n = resolve_parpool_workers()
%
%   Precedence:
%     1. PARPOOL_WORKERS          (explicit override)
%     2. SLURM_CPUS_PER_TASK      (sbatch --cpus-per-task)
%     3. SLURM_CPUS_ON_NODE       (--exclusive whole node)
%     4. SLURM_JOB_CPUS_PER_NODE
%     5. feature('numcores') / maxNumCompThreads('maximum')

  env_vars = {'PARPOOL_WORKERS', 'SLURM_CPUS_PER_TASK', ...
    'SLURM_CPUS_ON_NODE', 'SLURM_JOB_CPUS_PER_NODE'};
  for k = 1:numel(env_vars)
    nv = local_env_positive_int(env_vars{k});
    if ~isempty(nv)
      n = nv;
      return;
    end
  end

  try
    n = round(double(feature('numcores')));
  catch
    n = round(double(maxNumCompThreads('maximum')));
  end

  if ~isfinite(n) || n < 1
    n = 1;
  end
end

function n = local_env_positive_int(var_name)
  n = [];
  v = strtrim(getenv(var_name));
  if isempty(v)
    return;
  end
  nv = round(str2double(v));
  if isfinite(nv) && nv > 0
    n = nv;
  end
end
