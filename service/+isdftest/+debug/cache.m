function varargout = cache(cmd, varargin)
%CACHE  Internal persistent store for ISDF debug flags (set by init_from_config).
%
%   cache('init', config) 鈥?snapshot config.ISDF debug_* fields.
%   cache('get')        鈥?struct: initialized, checks, level, tags (cell).
%   cache('clear')      鈥?reset to uninitialized / checks off.

  persistent initialized checks level tags

  cmd = lower(strtrim(char(string(cmd))));

  switch cmd
    case 'init'
      cfg = [];
      if nargin >= 2
        cfg = varargin{1};
      end
      initialized = true;
      checks = false;
      level = 'error';
      tags = {};

      if isa(cfg, 'struct') && isfield(cfg, 'ISDF')
        s = cfg.ISDF;
        if isfield(s, 'debug_checks') && ~isempty(s.debug_checks)
          checks = logical(s.debug_checks);
        end
        if isfield(s, 'debug_level') && ~isempty(s.debug_level)
          v = lower(strtrim(char(string(s.debug_level))));
          if strcmp(v, 'warn') || strcmp(v, 'warning')
            level = 'warn';
          elseif strcmp(v, 'error')
            level = 'error';
          else
            warning('isdf:debug:badLevel', ...
              'Unknown config.isdftest.debug_level ''%s''; using ''error''.', char(string(s.debug_level)));
          end
        end
        if isfield(s, 'debug_tags') && ~isempty(s.debug_tags)
          tags = local_tags_as_cellstr(s.debug_tags);
        end
      end

    case 'clear'
      initialized = false;
      checks = false;
      level = 'error';
      tags = {};

    case 'get'
      if isempty(initialized) || ~logical(initialized)
        varargout{1} = struct('initialized', false, 'checks', false, ...
          'level', 'error', 'tags', {{}});
      else
        varargout{1} = struct('initialized', true, 'checks', logical(checks), ...
          'level', char(string(level)), 'tags', {tags});
      end

    otherwise
      error('isdf:debug:cache', 'Unknown cache command ''%s''.', cmd);
  end
end

function c = local_tags_as_cellstr(tags)
  if isstring(tags)
    c = cellstr(tags(:));
  elseif iscell(tags)
    c = tags(:);
  elseif ischar(tags)
    c = {char(tags)};
  else
    c = cellstr(string(tags));
  end
  keep = false(numel(c), 1);
  for i = 1:numel(c)
    keep(i) = strlength(strtrim(char(string(c{i})))) > 0;
  end
  c = c(keep);
  for i = 1:numel(c)
    c{i} = lower(strtrim(char(string(c{i}))));
  end
end
