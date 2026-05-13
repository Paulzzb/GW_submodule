function build_prod_mex()
%BUILD_PROD_MEX Build MEX backend for isdf.prod.
%
% Run from MATLAB:
%   isdf.build_prod_mex

  outdir = fileparts(mfilename('fullpath'));
  src = fullfile(outdir, 'prod_mex.c');
  mex('-R2018a', '-O', '-outdir', outdir, src, '-lmwblas');
end
