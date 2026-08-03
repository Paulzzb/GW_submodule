% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06 ZZ

function tf = on(varargin)
%ON  True when inline ISDF debug checks should run for this tag.
%
%   tf = isdf.debug.on(tag)
%       Uses flags snapshotted by isdf.debug.init_from_config(config) (e.g. from
%       service_driver). No config argument.
%
%   tf = isdf.debug.on(config, tag)   % optional overload
%       Reads config.ISDF directly (e.g. one-off scripts without init). Does not
%       update the session cache.
%
%   config.ISDF fields: debug_checks, debug_tags ([] = all tags when on),
%   debug_level is only used by react/level, not by on.

  tf = false;

  if nargin == 2 && isa(varargin{1}, 'struct') && isfield(varargin{1}, 'ISDF')
    tf = local_on_from_isdf(varargin{1}.ISDF, varargin{2});
    return
  end

  if nargin ~= 1
    error('isdf:debug:onArgs', 'Use isdf.debug.on(tag) or isdf.debug.on(config, tag).');
  end

  tag = varargin{1};
  if (isstring(tag) && ~isscalar(tag)) || isempty(tag)
    return
  end

  tag_norm = local_normalize_tag(tag);
  if strlength(tag_norm) == 0
    return
  end

  st = isdf.debug.cache('get');
  if ~st.initialized || ~st.checks
    return
  end

  tf = local_match_tags(st.tags, tag_norm);
end

function tf = local_on_from_isdf(s, tag)
  tf = false;
  tag_norm = local_normalize_tag(tag);
  if strlength(tag_norm) == 0
    return
  end
  if ~isfield(s, 'debug_checks') || isempty(s.debug_checks) || ~logical(s.debug_checks)
    return
  end
  tags = {};
  if isfield(s, 'debug_tags') && ~isempty(s.debug_tags)
    tags = local_tags_list(s.debug_tags);
  end
  tf = local_match_tags(tags, tag_norm);
end

function tf = local_match_tags(tags, tag_norm)
  if isempty(tags)
    tf = true;
    return
  end
  tf = false;
  for k = 1:numel(tags)
    if strcmpi(tag_norm, local_normalize_tag(tags{k}))
      tf = true;
      return
    end
  end
end

function t = local_normalize_tag(x)
  t = lower(strtrim(char(string(x))));
end

function c = local_tags_list(tags)
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
end
