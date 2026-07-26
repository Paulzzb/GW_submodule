function build_schur_rank1_mex()
%BUILD_SCHUR_RANK1_MEX Build rank-1 Schur update MEX kernel (double).
%
% Run this from MATLAB:
%   isdf.adaptive_double.build_schur_rank1_mex

  outdir = fileparts(mfilename('fullpath'));
  src = fullfile(outdir, 'isdf_schur_rank1_mex.c');
  local_clear_mex_stale(outdir, 'isdf_schur_rank1_mex');
  mex('-R2018a', '-O', '-outdir', outdir, src);
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
