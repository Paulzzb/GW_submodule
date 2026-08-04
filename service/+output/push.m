% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function cleanup = push(tag)
%PUSH  Push a module tag onto the display stack (for [tag] prefixes).
%
%   cleanup = output.push('qp_driver');
%   % ... tag auto-pops when cleanup goes out of scope

  if nargin < 1 || ~(ischar(tag) || isstring(tag))
    error('output:push:BadTag', 'tag must be a string.');
  end

  s = state_('get');
  s.module_stack{end+1} = char(string(tag)); %#ok<AGROW>
  state_('set', s);
  cleanup = onCleanup(@() local_pop());
end

function local_pop()
  s = state_('get');
  if ~isempty(s.module_stack)
    s.module_stack(end) = [];
    state_('set', s);
  end
end
