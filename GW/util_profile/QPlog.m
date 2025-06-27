function cleanup = QPlog(msg, level, setflag)
% QPlog Unified logging function for GW toolbox
%   QPlog(MSG, LEVEL)
%   LEVEL = 0: always shown
%   LEVEL = 1: normal progress messages
%   LEVEL = 2: verbose/debug messages

  persistent VERBOSE MODULE_STACK SHOWTAG ISLOGFILE LOGFILE_STATE

  cleanup = [];

  if isempty(VERBOSE); VERBOSE = 1; end % Default verbosity level
  if isempty(MODULE_STACK); MODULE_STACK = ''; end
  if isempty(SHOWTAG); SHOWTAG = 1; end
  if isempty(ISLOGFILE); ISLOGFILE = false; end

  if ischar(msg) || isstring(msg)
    msg = strtrim(msg);
  end



  % ----------------------------
  % Handle special calls
  % ----------------------------

  if nargin == 3
    switch setflag
      case 'verbose'
        % Setting verbosity level
        if isnumeric(msg)
          VERBOSE = msg;
        elseif ischar(msg) || isstring(msg)
          error('[QPlog] VERBOSE expects to be integer.');
        end
        return
      case 'showtag'
        if islogical(msg)
          SHOWTAG = msg;
        else
          error('[QPlog] showtag expects logical true/false.');
        end
        return
      case 'push'
        if ischar(msg) || isstring(msg)
          MODULE_STACK{end+1} = msg;
          cleanup = onCleanup(@() QPlog_pop());
        else
          error('[QPlog] MODULE name expects to be string.');
        end
        return     
      case 'logfile'
        if ischar(msg) || isstring(msg)
          ISLOGFILE = true;
          [fid, errmsg] = fopen(msg, 'a');
          if fid < 0
            error(['Fail to open logfile: %s', errmsg], logfile);
          end
          LOGFILE_STATE.HANDLE = fid;
          cleanupObj = onCleanup(@()fclose(LOGFILE_STATE.HANDLE));
          LOGFILE_STATE.CLEANUP = cleanupObj;
        else
          error('[QPlog] MODULE name expects to be string.');
        end
        return
      otherwise
        error('[QPlog] unknown flag "%s".', setflag);
    end
  end

  % ----------------------------
  % Default verbosity level if not provided
  % ----------------------------
  if (nargin < 2)
    level = 1;
  end
  

  % ----------------------------
  % Conditional output based on verbosity
  % ----------------------------
  if level <= VERBOSE
    if SHOWTAG && ~isempty(MODULE_STACK)
      prefix = sprintf('[%s] ', MODULE_STACK{end});
    else
      prefix = '';
    end
    if ISLOGFILE
      fprintf(LOGFILE_STATE.HANDLE, '%s%s\n', prefix, msg);
    else
      fprintf('%s%s\n', prefix, msg);
    end
  end

  function QPlog_pop()
    if ~isempty(MODULE_STACK)
      MODULE_STACK(end) = [];
    end
  end
end

