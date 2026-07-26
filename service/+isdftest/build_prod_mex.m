function build_prod_mex()
%BUILD_PROD_MEX Build MEX backend for isdftest.prod (double).
%
% Run from MATLAB:
%   isdftest.build_prod_mex

  outdir = fileparts(mfilename('fullpath'));
  src = fullfile(outdir, 'prod_mex.c');
  local_clear_mex_stale(outdir, 'prod_mex');
  mex('-R2018a', '-O', '-outdir', outdir, src, '-lmwblas');
end

function local_clear_mex_stale(outdir, stem)
  stale = dir(fullfile(outdir, [stem '.mex*']));
  for k = 1:numel(stale)
    f = fullfile(outdir, stale(k).name);
    if isfile(f)
      delete(f);
    end
  end
end
