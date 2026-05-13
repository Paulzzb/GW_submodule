function build_schur_rank1_mex()
%BUILD_SCHUR_RANK1_MEX Build rank-1 Schur update MEX kernel.
%
% Run this from MATLAB:
%   isdf.adaptive.build_schur_rank1_mex

  outdir = fileparts(mfilename('fullpath'));
  src = fullfile(outdir, 'isdf_schur_rank1_mex.c');
  mex('-R2018a', '-O', '-outdir', outdir, src);
end
