% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% Persistent pool of isdf_m (N_MAX = 10), Yambo FFT(FFT_N_max) style.
% Commands:
%   'get'           [, id]   -> isdf_m at id or at current id
%   'save2mod', obj [, id]  -> write obj (default: current id)
%   'add', desc             -> first free cell, placeholder, set current id, return id
%   'select', id            -> set current id
%   'current'               -> current id (int32)
%   'nmax'                  -> N_MAX (int32)
%   'list'                  -> struct array id/desc/assigned/empty
%   'free'                  -> clear entire pool
%   'free', id              -> clear one id; adjust current if needed

function out = manager(cmd, varargin)
  persistent pool current_id pool_init

  N_MAX = int32(10);

  if isempty(pool_init)
    pool = cell(1, double(N_MAX));
    current_id = int32(1);
    pool_init = true;
  end

  if nargin == 0
    cmd = 'get';
  end

  switch lower(cmd)
    case 'get'
      pid = current_id;
      if nargin >= 2 && ~isempty(varargin{1})
        pid = int32(varargin{1});
      end
      isdf_manager_validate_id(pid, N_MAX);
      if isempty(pool{double(pid)})
        error('ISDF:manager:EmptyId', 'ISDF id %d is empty.', pid);
      end
      out = pool{double(pid)};
      return

    case 'save2mod'
      if isempty(varargin)
        error('ISDF:manager:Save2modArgs', 'ISDF::save2mod requires data.');
      end
      candidate = varargin{1};
      if ~(isa(candidate, 'isdftest.base.isdftest_m'))
        error('ISDF:manager:Save2modType', 'ISDF::save2mod requires isdftest.base.isdftest_m.');
      end
      pid = current_id;
      if nargin >= 3 && ~isempty(varargin{2})
        pid = int32(varargin{2});
      end
      isdf_manager_validate_id(pid, N_MAX);
      candidate.id = pid;
      pool{double(pid)} = candidate;
      out = int32(0);
      return

    case 'add'
      desc_str = "";
      if nargin >= 2 && ~isempty(varargin{1})
        desc_str = string(varargin{1});
      end
      for i = 1:double(N_MAX)
        if isempty(pool{i})
          obj = isdftest.base.isdftest_m();
          obj.desc = desc_str;
          obj.id = int32(i);
          obj.allocated = true;
          obj.assigned = false;
          pool{i} = obj;
          current_id = int32(i);
          out = int32(i);
          return;
        end
      end
      error('ISDF:manager:PoolFull', 'ISDF pool full (N_MAX = %d).', N_MAX);

    case 'select'
      if isempty(varargin)
        error('ISDF:manager:SelectArgs', 'ISDF::select requires id.');
      end
      pid = int32(varargin{1});
      isdf_manager_validate_id(pid, N_MAX);
      if isempty(pool{double(pid)})
        error('ISDF:manager:SelectEmpty', 'ISDF id %d is empty.', pid);
      end
      current_id = pid;
      out = int32(0);
      return

    case 'current'
      out = current_id;
      return

    case 'nmax'
      out = N_MAX;
      return

    case 'list'
      L = struct('id', cell(1, double(N_MAX)), 'desc', [], 'assigned', [], 'empty', []);
      for i = 1:double(N_MAX)
        L(i).id = int32(i);
        if isempty(pool{i})
          L(i).desc = "";
          L(i).assigned = false;
          L(i).empty = true;
        else
          L(i).desc = pool{i}.desc;
          L(i).assigned = pool{i}.assigned;
          L(i).empty = false;
        end
      end
      out = L;
      return

    case 'free'
      if nargin < 2 || isempty(varargin{1})
        pool = cell(1, double(N_MAX));
        current_id = int32(1);
        out = int32(0);
        return;
      end
      pid = int32(varargin{1});
      isdf_manager_validate_id(pid, N_MAX);
      pool{double(pid)} = [];
      if current_id == pid
        current_id = int32(1);
        for j = 1:double(N_MAX)
          if ~isempty(pool{j})
            current_id = int32(j);
            break;
          end
        end
      end
      out = int32(0);
      return

    otherwise
      error('ISDF:manager:UnknownCommand', 'Unknown command: %s', cmd);
  end
end

function isdf_manager_validate_id(pid, N_MAX)
  if ~(isnumeric(pid) && isscalar(pid) && pid == floor(double(pid)))
    error('ISDF:manager:BadId', 'id must be a scalar integer.');
  end
  pid = double(pid);
  if pid < 1 || pid > double(N_MAX)
    error('ISDF:manager:BadId', 'id must be in 1..%d.', double(N_MAX));
  end
end
