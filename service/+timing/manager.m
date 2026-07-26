% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function out = manager(cmd, varargin)
  persistent data

  if nargin == 0
    cmd = 'get';
  end

  switch lower(cmd)
    case 'get'
      if isempty(data)
        error('timing not initialized. Call timing.driver first.');
      end
      out = data;
      return

    case 'save2mod'
      if isempty(varargin)
        error('timing::save2mod requires one input argument.');
      end

      candidate = varargin{1};
      if ~(isa(candidate, 'timing.base.timing_m'))
        error('timing::save2mod requires input as a timing.base.timing_m object.');
      end

      data = candidate;
      out = int32(0);

    case 'free'
      data = [];
      out = int32(0);

    otherwise
      error('timing::manager:UnknownCommand', 'Unknown command: %s', cmd);
  end
end
