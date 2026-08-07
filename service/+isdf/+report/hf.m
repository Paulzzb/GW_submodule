% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ

function fpath = hf(report)
%HF  Write HF validation text report (from isdf_validate_HF).
%
%   fpath = isdf.report.hf(report)
%
% Always under filename_map().isdf_report_dir / sprintf(hf_report, id).
% Expects plain fields from isdf_validate_HF (no MATLAB table).
% See +report/NAMING.md.

  def = filename_map();
  id = double(report.isdf_id);
  of = sprintf(def.hf_report, id);
  how = ['o ' of];
  fpath = fullfile(def.isdf_report_dir, of);
  output.open(of, fpath, 'w');

  output.msg(how, '=== ISDF validate HF report (isdf id = %d) ===', id);
  output.msg(how, 'Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
  if isfield(report, 'desc') && strlength(string(report.desc)) > 0
    output.msg(how, 'description: %s', char(string(report.desc)));
  end
  if isfield(report, 'interp_scheme')
    output.msg(how, 'interp_scheme: %s', char(string(report.interp_scheme)));
  end
  if isfield(report, 'nisdf')
    output.msg(how, 'nisdf: %d', double(report.nisdf));
  end
  output.msg(how, '');

  output.msg(how, '--- Band-resolved summary (E_HF vs E_HF_ISDF) ---');
  output.msg(how, 'Definition: per (ib, ik_ibz, ispin), E_HF and E_HF_ISDF sum over all ob paths.');
  output.msg(how, '');
  for k = 1:numel(report.band_items)
    output.msg(how, '  %-48s  %s', report.band_items{k}, report.band_vals{k});
  end

  output.msg(how, '');
  output.msg(how, '--- Preview (first min(N, preview_rows) ob samples; stats use full online accumulation) ---');
  pr = report.preview;
  nprev = 0;
  if isstruct(pr) && isfield(pr, 'Ex_t')
    nprev = numel(pr.Ex_t);
  end
  if nprev > 0
    output.msg(how, 'preview_rows budget: %d', double(report.preview_rows));
    output.msg(how, '');
    output.msg(how, '  %12s %12s %12s %12s %18s %18s', ...
      'Ex_t', 'Ex_ISDF', 'Diff', 'AbsDiff', 'Diff/|Ex_t|', 'AbsDiff/|Ex_t|');
    for i = 1:nprev
      output.msg(how, '  %12.6g %12.6g %12.6g %12.6g %18.6g %18.6g', ...
        pr.Ex_t(i), pr.Ex_ISDF(i), pr.Diff(i), pr.AbsDiff(i), ...
        pr.RelSigned(i), pr.RelAbs(i));
    end
  else
    output.msg(how, '(no ob samples: all skipped by occupation threshold)');
  end

  output.msg(how, '');
  output.msg(how, '--- Statistics (Mean / Std / Var / Max), ob layer ---');
  output.msg(how, 'N (total ob samples) = %d', double(report.n_samples));
  output.msg(how, 'N_rel (samples with |Ex_t| >= tol for relative stats) = %d', ...
    double(report.n_samples_rel));
  output.msg(how, '');
  vars = report.stats_vars;
  rows = report.stats_rows;
  M = report.stats_M;
  hdr = sprintf('  %-6s', '');
  for j = 1:numel(vars)
    hdr = [hdr, sprintf(' %14s', vars{j})]; %#ok<AGROW>
  end
  output.msg(how, '%s', hdr);
  for i = 1:numel(rows)
    line = sprintf('  %-6s', rows{i});
    for j = 1:size(M, 2)
      line = [line, sprintf(' %14.6g', M(i, j))]; %#ok<AGROW>
    end
    output.msg(how, '%s', line);
  end
  output.msg(how, '');
  output.msg(how, 'Note: Diff = Ex_t - Ex_ISDF; relative columns use only samples with |Ex_t| >= tol (N_rel may differ from N).');

  output.msg(how, '');
  output.msg(how, '--- Global energy sums (accumulated over ob paths) ---');
  output.msg(how, 'Sum Ex_t^2 (direct):      %.12e', report.Esum2);
  output.msg(how, 'Sum Ex_ISDF^2:            %.12e', report.EsumISDF2);
  output.msg(how, 'Sum (Ex_t - Ex_ISDF)^2:   %.12e', report.DiffEsum2);
  output.msg(how, '');
  output.msg(how, '=== end of report ===');

  output.close(of);
end
