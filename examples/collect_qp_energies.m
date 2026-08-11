% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/08 ZZ

function T = collect_qp_energies(cases_root, case_names)
%COLLECT_QP_ENERGIES  Gather Eqp0 (and related columns) from case qp*.dat files.
%
%   T = collect_qp_energies()
%   T = collect_qp_energies(cases_root)
%   T = collect_qp_energies(cases_root, case_names)
%
%   Scans examples/cases/<name>/qp.dat (and qp_*.dat) written by qp.fout.
%   Supports:
%     frequency_dependence == 2   — 2 lines/band; uses Re Eqp0
%     frequency_dependence == -2  — 1 line/band; uses Eqp0
%
%   Returns a table with columns:
%     case, file, band, Emf, Eo, X, SXx, CH, Sig, Vxc, Eqp0
%   and prints a compact Eqp0 comparison (bands x cases) to the command window.

  if nargin < 1 || isempty(cases_root)
    here = fileparts(mfilename('fullpath'));
    cases_root = fullfile(here, 'cases');
  end
  if nargin < 2 || isempty(case_names)
    case_names = { ...
      'gamma_ff_isdf', ...
      'gamma_ff_dir', ...
      'gamma_cohsex_isdf', ...
      'gamma_cohsex_dir' ...
      };
  end

  rows = {};
  for ic = 1:numel(case_names)
    cname = case_names{ic};
    cdir = fullfile(cases_root, cname);
    if ~isfolder(cdir)
      fprintf('SKIP (missing dir): %s\n', cname);
      continue
    end
    hits = [dir(fullfile(cdir, 'qp.dat')); dir(fullfile(cdir, 'qp_*.dat'))];
    if isempty(hits)
      fprintf('SKIP (no qp*.dat): %s\n', cname);
      continue
    end
    for ih = 1:numel(hits)
      if hits(ih).isdir
        continue
      end
      fpath = fullfile(cdir, hits(ih).name);
      bands = local_parse_qp_dat(fpath);
      if isempty(bands)
        fprintf('SKIP (parse empty): %s\n', fpath);
        continue
      end
      for ib = 1:numel(bands)
        b = bands(ib);
        rows(end + 1, :) = {cname, hits(ih).name, b.n, ...
          b.Emf, b.Eo, b.X, b.SXx, b.CH, b.Sig, b.Vxc, b.Eqp0}; %#ok<AGROW>
      end
      fprintf('OK  %-28s  %s  (%d bands)\n', cname, hits(ih).name, numel(bands));
    end
  end

  if isempty(rows)
    T = table();
    warning('collect_qp_energies:NoData', 'No qp*.dat tables found.');
    return
  end

  T = cell2table(rows, 'VariableNames', ...
    {'case', 'file', 'band', 'Emf', 'Eo', 'X', 'SXx', 'CH', 'Sig', 'Vxc', 'Eqp0'});

  % Compact Eqp0 pivot: rows=band, columns=case (qp.dat only when several files).
  Tqp = T(strcmp(T.file, 'qp.dat'), :);
  if isempty(Tqp)
    Tqp = T;
  end
  bands = unique(Tqp.band, 'stable');
  cases_u = unique(Tqp.case, 'stable');
  M = nan(numel(bands), numel(cases_u));
  for i = 1:height(Tqp)
    ib = find(bands == Tqp.band(i), 1);
    ic = find(strcmp(cases_u, Tqp.case{i}), 1);
    M(ib, ic) = Tqp.Eqp0(i);
  end

  fprintf('\n=== Re Eqp0 / Eqp0 (eV) ===\n');
  hdr = sprintf('%6s', 'band');
  for ic = 1:numel(cases_u)
    hdr = [hdr, sprintf('  %16s', cases_u{ic})]; %#ok<AGROW>
  end
  fprintf('%s\n', hdr);
  for ib = 1:numel(bands)
    line = sprintf('%6d', bands(ib));
    for ic = 1:numel(cases_u)
      if isnan(M(ib, ic))
        line = [line, sprintf('  %16s', '—')]; %#ok<AGROW>
      else
        line = [line, sprintf('  %16.6f', M(ib, ic))]; %#ok<AGROW>
      end
    end
    fprintf('%s\n', line);
  end

  out_csv = fullfile(cases_root, '..', 'qp_energies_collect.csv');
  try
    writetable(T, out_csv);
    fprintf('\nWrote %s\n', out_csv);
  catch ME
    warning('collect_qp_energies:CsvFailed', '%s', ME.message);
  end
end

function bands = local_parse_qp_dat(fpath)
  bands = struct('n', {}, 'Emf', {}, 'Eo', {}, 'X', {}, 'SXx', {}, ...
    'CH', {}, 'Sig', {}, 'Vxc', {}, 'Eqp0', {});
  fid = fopen(fpath, 'r');
  if fid < 0
    return
  end
  cleaner = onCleanup(@() fclose(fid));
  lines = {};
  while true
    ln = fgetl(fid);
    if ~ischar(ln)
      break
    end
    lines{end + 1} = ln; %#ok<AGROW>
  end
  if numel(lines) < 2
    return
  end

  is_ff = contains(lines{1}, 'Re Eqp0') || contains(lines{1}, 'Re SX-X');
  i = 1;
  % skip header lines
  while i <= numel(lines) && (contains(lines{i}, 'Emf') || contains(lines{i}, 'Im SX') ...
      || contains(lines{i}, 'SX-X') && ~local_starts_with_int(lines{i}))
    i = i + 1;
  end

  while i <= numel(lines)
    ln = strtrim(lines{i});
    if isempty(ln)
      i = i + 1;
      continue
    end
    vals = sscanf(ln, '%f');
    if numel(vals) < 8
      i = i + 1;
      continue
    end
    b = struct();
    b.n = vals(1);
    b.Emf = vals(2);
    b.Eo = vals(3);
    b.X = vals(4);
    b.SXx = vals(5);
    b.CH = vals(6);
    b.Sig = vals(7);
    b.Vxc = vals(8);
    b.Eqp0 = vals(9);
    bands(end + 1) = b; %#ok<AGROW>
    i = i + 1;
    if is_ff && i <= numel(lines)
      % skip imaginary continuation line
      ln2 = strtrim(lines{i});
      v2 = sscanf(ln2, '%f');
      if ~isempty(v2) && (numel(v2) < 8 || ~local_starts_with_int(ln2))
        i = i + 1;
      end
    end
  end
end

function tf = local_starts_with_int(ln)
  tf = ~isempty(regexp(strtrim(ln), '^\d+', 'once'));
end
