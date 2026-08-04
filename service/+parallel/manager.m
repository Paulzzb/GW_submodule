% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06

function out = manager(cmd, varargin)
  persistent data

  if nargin == 0
    cmd = 'get';
  end

  switch lower(cmd)
    case 'get'
      if isempty(data)
        error('parallel not initialized. Call parallel.driver first.');
      end
      out = data;
      return

    case 'save2mod'
      if isempty(varargin)
        error('parallel::save2mod requires one input argument.');
      end

      candidate = varargin{1};
      if ~(isa(candidate, 'parallel.base.parallel_m'))
        error('parallel::save2mod requires input as a parallel.base.parallel_m object.');
      end

      data = candidate;
      out = int32(0);

    case 'free'
      data = [];
      out = int32(0);

    otherwise
      error('parallel::manager:UnknownCommand', 'Unknown command: %s', cmd);
  end
end
