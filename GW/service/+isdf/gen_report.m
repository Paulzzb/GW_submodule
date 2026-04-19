% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/16

function fpath = gen_report(report, outDir)
% isdf.gen_report  Write HF validation text report produced by isdf_validate_HF to a file.
%
%   fpath = isdf.gen_report(report)
%   fpath = isdf.gen_report(report, outDir)   % optional directory (default: pwd)
%
% Required fields on report:
%   isdf_id          — ISDF pool id (integer)
%   band_summary     — table with columns Item, Result
%   preview          — table (may be empty)
%   stats            — table of ob-layer statistics
%   n_samples        — total ob sample count N
%   n_samples_rel    — count used for relative-error streaming stats
%   preview_rows     — preview row budget (min(N, preview_rows))
%   Esum2, EsumISDF2, DiffEsum2 — accumulated sum of squares (same as validate_HF returns)
%
% Optional: desc, interp_scheme, nisdf
%
% Output file:
%   <outDir>/isdf_validate_HF_id<isdf_id>.txt   (outDir defaults to pwd)

  if ~isstruct(report) || ~isfield(report, 'isdf_id')
    error('isdf:gen_report:BadReport', 'report must be a struct with field isdf_id.');
  end

  if nargin < 2
    outDir = [];
  end
  if isempty(outDir)
    outDir = pwd;
  else
    outDir = char(string(outDir));
  end
  if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
  end

  id = double(report.isdf_id);
  fname = sprintf('isdf_validate_HF_id%d.txt', id);
  fpath = fullfile(outDir, fname);

  fid = fopen(fpath, 'w');
  if fid < 0
    error('isdf:gen_report:OpenFailed', 'Cannot open file for write: %s', fpath);
  end
  oc = onCleanup(@() fclose(fid));

  fprintf(fid, '=== ISDF validate HF report (isdf id = %d) ===\n', id);
  fprintf(fid, 'Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
  if isfield(report, 'desc') && strlength(string(report.desc)) > 0
    fprintf(fid, 'description: %s\n', char(string(report.desc)));
  end
  if isfield(report, 'interp_scheme')
    fprintf(fid, 'interp_scheme: %s\n', char(string(report.interp_scheme)));
  end
  if isfield(report, 'nisdf')
    fprintf(fid, 'nisdf: %d\n', double(report.nisdf));
  end
  fprintf(fid, '\n');

  fprintf(fid, '--- Band-resolved summary (E_HF vs E_HF_ISDF) ---\n');
  fprintf(fid, 'Definition: per (ib, ik_ibz, ispin), E_HF and E_HF_ISDF sum over all ob paths.\n\n');
  fprintf(fid, '%s\n', strtrim(isdf_gen_report_table2char(report.band_summary)));

  fprintf(fid, '\n--- Preview (first min(N, preview_rows) ob samples; stats use full online accumulation) ---\n');
  pr = report.preview;
  if istable(pr) && height(pr) > 0
    fprintf(fid, 'preview_rows budget: %d\n\n', double(report.preview_rows));
    fprintf(fid, '%s\n', strtrim(isdf_gen_report_table2char(pr)));
  else
    fprintf(fid, '(no ob samples: all skipped by occupation threshold)\n');
  end

  fprintf(fid, '\n--- Statistics (Mean / Std / Var / Max), ob layer ---\n');
  fprintf(fid, 'N (total ob samples) = %d\n', double(report.n_samples));
  fprintf(fid, 'N_rel (samples with |Ex_t| >= tol for relative stats) = %d\n\n', double(report.n_samples_rel));
  fprintf(fid, '%s\n', strtrim(isdf_gen_report_table2char(report.stats)));
  fprintf(fid, '\nNote: Diff = Ex_t - Ex_ISDF; relative columns use only samples with |Ex_t| >= tol (N_rel may differ from N).\n');

  fprintf(fid, '\n--- Global energy sums (accumulated over ob paths) ---\n');
  fprintf(fid, 'Sum Ex_t^2 (direct):      %.12e\n', report.Esum2);
  fprintf(fid, 'Sum Ex_ISDF^2:            %.12e\n', report.EsumISDF2);
  fprintf(fid, 'Sum (Ex_t - Ex_ISDF)^2:   %.12e\n', report.DiffEsum2);
  fprintf(fid, '\n=== end of report ===\n');

  % onCleanup closes fid
end

function s = isdf_gen_report_table2char(tbl)
  if ~istable(tbl) || isempty(tbl) || width(tbl) == 0
    s = '';
    return
  end
  s = evalc('disp(tbl)');
end
