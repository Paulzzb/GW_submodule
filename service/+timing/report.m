% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Print a fixed-width timing report (global + internal clocks).
% Layout: horizontal rule / title + column headers / rule / data rows + TOTAL / rule.
% If logfile is non-empty (argument or tm.live.report_logfile), append to that file;
% otherwise print to the command window (fid 1).
%
%   timing.report()
%   timing.report('logfile', '/path/to/run.log')
%   timing.report('logfile', '', 'tm', tm_obj)   % use given timing_m, no save
%
% MATLAB R2008a: addParamValue, char paths only.

function report(varargin)
  p = inputParser;
  p.FunctionName = mfilename;
  p.addParamValue('logfile', '', @ischar);
  p.addParamValue('tm', [], @(x) isempty(x) || isa(x, 'timing.base.timing_m'));
  p.parse(varargin{:});

  logfile = strtrim(p.Results.logfile);
  tm = p.Results.tm;

  if isempty(tm)
    tm = timing.get();
  end

  if isempty(logfile)
    logfile = strtrim(char(tm.live.report_logfile));
  end

  if isempty(logfile)
    fid = 1;
    must_close = false;
  else
    fid = fopen(logfile, 'a');
    if fid < 0
      error('timing:report:OpenFailed', 'Could not open log file for append: %s', logfile);
    end
    must_close = true;
  end

  try
    lines = timing_report_lines(tm);
    for i = 1:numel(lines)
      fprintf(fid, '%s\n', lines{i});
    end
    if must_close
      fclose(fid);
    end
  catch ME
    if must_close && fid > 0
      fclose(fid);
    end
    rethrow(ME);
  end
end

function lines = timing_report_lines(tm)
  sep = repmat('-', 1, 80);

  rows = timing_report_collect_rows(tm);
  nw = timing_report_name_width(rows);

  title_line = sprintf('%-*s %8s %14s %8s', nw, 'Clock (list/name)', 'Calls', 'CPU (s)', 'Run');

  % Three-rule layout: ---- / titles / ---- / info (rows + total) / ----
  lines = cell(0, 1);
  lines{end + 1, 1} = sep; %#ok<AGROW>
  lines{end + 1, 1} = 'GW TIMING REPORT'; %#ok<AGROW>
  lines{end + 1, 1} = title_line; %#ok<AGROW>
  lines{end + 1, 1} = sep; %#ok<AGROW>

  total_cpu = 0;
  for r = 1:numel(rows)
    lab = timing_report_fit_label(rows(r).label, nw);
    nc = rows(r).calls;
    tc = rows(r).cpu;
    rn = rows(r).run;
    total_cpu = total_cpu + tc;
    lines{end + 1, 1} = sprintf('%-*s %8d %14.6f %8s', nw, lab, nc, tc, rn); %#ok<AGROW>
  end

  if numel(rows) == 0
    lines{end + 1, 1} = '(no allocated clocks)'; %#ok<AGROW>
  end

  tot_lbl = timing_report_fit_label('TOTAL (listed)', nw);
  lines{end + 1, 1} = sprintf('%-*s %8s %14.6f %8s', nw, tot_lbl, '--', total_cpu, '--'); %#ok<AGROW>
  lines{end + 1, 1} = sep; %#ok<AGROW>
end

function rows = timing_report_collect_rows(tm)
  rows = [];
  rows = timing_report_append_list(rows, tm.internal_list);
  rows = timing_report_append_list(rows, tm.global_list);
end

function rows = timing_report_append_list(rows, list)
  if isempty(list) || ~list.alloc
    return;
  end
  lname = char(list.name);
  n = double(list.nclock);
  for i = 1:n
    c = list.clocks(i);
    if ~c.alloc
      continue;
    end
    nm = char(c.name);
    lab = [lname, '/', nm];
    r.label = lab;
    r.calls = double(c.call_number);
    r.cpu = double(c.total_time);
    if c.running
      r.run = 'yes';
    else
      r.run = 'no';
    end
    if isempty(rows)
      rows = r;
    else
      rows = [rows; r]; %#ok<AGROW>
    end
  end
end

function nw = timing_report_name_width(rows)
  nw = numel('Clock (list/name)');
  for i = 1:numel(rows)
    nw = max(nw, numel(rows(i).label));
  end
  nw = max(nw, numel('TOTAL (listed)'));
  if nw < 24
    nw = 24;
  end
  if nw > 56
    nw = 56;
  end
end

function out = timing_report_fit_label(lab, nw)
  lab = char(lab);
  if numel(lab) <= nw
    out = lab;
    return;
  end
  if nw <= 3
    out = lab(1:nw);
    return;
  end
  out = [lab(1:nw - 3), '...'];
end
