% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function fpath = hf(report, outDir)
%HF  Write HF validation text report (from isdf_validate_HF).
%
%   fpath = isdf.report.hf(report)
%   fpath = isdf.report.hf(report, outDir)   % default: pwd
%
% Disk name: filename_map().hf_report -> o-ISDF_HF_id%d
% Also writes legacy isdf_validate_HF_id%d.txt for one compatibility round.
% See +report/NAMING.md.

  if nargin < 2 || isempty(outDir)
    outDir = pwd;
  else
    outDir = char(string(outDir));
  end
  if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
  end

  id = double(report.isdf_id);
  def = filename_map();
  fpath = fullfile(outDir, sprintf(def.hf_report, id));
  output.open('hf_report', fpath, 'w');

  output.msg('o hf_report', '=== ISDF validate HF report (isdf id = %d) ===', id);
  output.msg('o hf_report', 'Generated: %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
  if isfield(report, 'desc') && strlength(string(report.desc)) > 0
    output.msg('o hf_report', 'description: %s', char(string(report.desc)));
  end
  if isfield(report, 'interp_scheme')
    output.msg('o hf_report', 'interp_scheme: %s', char(string(report.interp_scheme)));
  end
  if isfield(report, 'nisdf')
    output.msg('o hf_report', 'nisdf: %d', double(report.nisdf));
  end
  output.msg('o hf_report', '');

  output.msg('o hf_report', '--- Band-resolved summary (E_HF vs E_HF_ISDF) ---');
  output.msg('o hf_report', 'Definition: per (ib, ik_ibz, ispin), E_HF and E_HF_ISDF sum over all ob paths.');
  output.msg('o hf_report', '');
  local_msg_block(local_table2char(report.band_summary));

  output.msg('o hf_report', '');
  output.msg('o hf_report', '--- Preview (first min(N, preview_rows) ob samples; stats use full online accumulation) ---');
  pr = report.preview;
  if istable(pr) && height(pr) > 0
    output.msg('o hf_report', 'preview_rows budget: %d', double(report.preview_rows));
    output.msg('o hf_report', '');
    local_msg_block(local_table2char(pr));
  else
    output.msg('o hf_report', '(no ob samples: all skipped by occupation threshold)');
  end

  output.msg('o hf_report', '');
  output.msg('o hf_report', '--- Statistics (Mean / Std / Var / Max), ob layer ---');
  output.msg('o hf_report', 'N (total ob samples) = %d', double(report.n_samples));
  output.msg('o hf_report', 'N_rel (samples with |Ex_t| >= tol for relative stats) = %d', double(report.n_samples_rel));
  output.msg('o hf_report', '');
  local_msg_block(local_table2char(report.stats));
  output.msg('o hf_report', '');
  output.msg('o hf_report', 'Note: Diff = Ex_t - Ex_ISDF; relative columns use only samples with |Ex_t| >= tol (N_rel may differ from N).');

  output.msg('o hf_report', '');
  output.msg('o hf_report', '--- Global energy sums (accumulated over ob paths) ---');
  output.msg('o hf_report', 'Sum Ex_t^2 (direct):      %.12e', report.Esum2);
  output.msg('o hf_report', 'Sum Ex_ISDF^2:            %.12e', report.EsumISDF2);
  output.msg('o hf_report', 'Sum (Ex_t - Ex_ISDF)^2:   %.12e', report.DiffEsum2);
  output.msg('o hf_report', '');
  output.msg('o hf_report', '=== end of report ===');

  output.close('hf_report');

  % Compatibility round: keep old basename for existing tests / parsers.
  fpath_legacy = fullfile(outDir, sprintf('isdf_validate_HF_id%d.txt', id));
  copyfile(fpath, fpath_legacy);
end

function local_msg_block(s)
  if isempty(s)
    return
  end
  lines = regexp(char(string(s)), '\r\n|\n|\r', 'split');
  for k = 1:numel(lines)
    output.msg('o hf_report', '%s', lines{k});
  end
end

function s = local_table2char(tbl)
  if ~istable(tbl) || isempty(tbl) || width(tbl) == 0
    s = '';
    return
  end
  s = strtrim(evalc('disp(tbl)'));
end
