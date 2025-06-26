function QPlog(msg, level, setflag)
% GWlog Unified logging function for GW toolbox
%   QPlog(MSG, LEVEL)
%   LEVEL = 0: always shown
%   LEVEL = 1: normal progress messages
%   LEVEL = 2: verbose/debug messages

  persistent VERBOSE MODULE SHOWTAG

  if isempty(VERBOSE); VERBOSE = 1; end % Default verbosity level
  if isempty(MODULE); MODULE = ''; end
  if isempty(SHOWTAG); SHOWTAG = 1; end

  if ischar(msg) || isstring(msg)
    msg = strtrim(msg);
  end

  if isnumeric(msg)
    VERBOSE = msg;
    return;
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
          error('GWlog: VERBOSE expects to be integer.');
        end
        return
      case 'module'
        if ischar(msg) || isstring(msg)
          MODULE = char(msg);
        else
          error('GWlog: MODULE expects to be string.');
        end
        return
      case 'showtag'
        if islogical(msg)
          SHOWTAG = msg;
        else
          error('GWlog: showtag expects logical true/false.');
        end
        return
      otherwise
        error('GWlog: unknown flag "%s".', setflag);
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
    if SHOWTAG && ~isempty(MODULE)
      prefix = sprintf('[%s] ', MODULE);
    else
      prefix = '';
    end
    fprintf('%s%s\n', prefix, msg);
  end

end

