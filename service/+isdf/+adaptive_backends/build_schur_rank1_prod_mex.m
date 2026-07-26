function build_schur_rank1_prod_mex()
%BUILD_SCHUR_RANK1_PROD_MEX Build rank-1 Schur update MEX with fused prod.
%
% Run this from MATLAB:
%   isdf.adaptive.build_schur_rank1_prod_mex

  outdir = fileparts(mfilename('fullpath'));
  src = fullfile(outdir, 'isdf_schur_rank1_prod_mex.c');
  mex('-R2018a', '-O', '-outdir', outdir, src, '-lmwblas');
end
