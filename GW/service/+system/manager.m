% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/26

function out = manager(cmd, varargin)
  persistent system_data

  if nargin == 0
    cmd = 'get';
  end

  switch lower(cmd)
    case 'get'
      if isempty(system_data)
        error('system not initialized. Call system.save2mod or system.driver first.');
      end
      out = system_data;
      return

    case 'save2mod'
      if isempty(varargin)
        error('system::save2mod requires one input argument.');
      end

      candidate = varargin{1};
      if ~(isa(candidate, 'system.base.system_m'))
        error('system::save2mod requires input as a system.base.system_m object.');
      end

      system_data = candidate;
      out = int32(0);

    case 'free'
      system_data = [];
      out = int32(0);

    otherwise
      error('Unknown command: %s', cmd);
  end
end