% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function emit_(flags, line)
%EMIT_  Write one (already formatted) line to the destinations in flags.

  s = state_('get');

  if flags.verbose > s.verbose
    return
  end

  if flags.blank_before
    local_write_dest(s, flags, '');
  end

  local_write_dest(s, flags, line);

  if flags.blank_after
    local_write_dest(s, flags, '');
  end
end

function local_write_dest(s, flags, line)
  % Named output file
  if ~isempty(flags.of_name)
    fid = local_named_fid(s, flags.of_name);
    if fid > 0
      fprintf(fid, '%s\n', line);
    end
    return
  end

  want_screen = flags.screen && s.write_to_screen;
  want_report = flags.report && s.write_to_report && s.report_fid > 0;
  want_log = flags.log && s.write_to_log && s.log_fid > 0;

  % Legacy set_logfile mode: screen-bound messages go exclusively to the log file.
  if s.legacy_logfile_only && want_screen && ~want_log
    want_screen = false;
    if s.log_fid > 0
      want_log = true;
    else
      want_screen = true;
    end
  end

  if want_screen
    fprintf('%s\n', line);
  end
  if want_report
    fprintf(s.report_fid, '%s\n', line);
  end
  if want_log
    fprintf(s.log_fid, '%s\n', line);
  end
end

function fid = local_named_fid(s, name)
  fid = -1;
  for i = 1:numel(s.named)
    if strcmp(s.named(i).name, name)
      fid = s.named(i).fid;
      return
    end
  end
end
