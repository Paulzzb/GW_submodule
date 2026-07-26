% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/31

function out = manager(cmd, varargin)
  % Persistent cache for the pair-symmetry model.
  % Mirrors the manager pattern used by other service modules.
  persistent pair_symm_data

  if nargin == 0
    cmd = 'get';
  end

  switch lower(cmd)
    case 'get'
      if isempty(pair_symm_data)
        error('pair_symmetry not initialized. Call pair_symmetry.save2mod or pair_symmetry.driver first.');
      end
      out = pair_symm_data;
      return

    case 'save2mod'
      if isempty(varargin)
        error('pair_symmetry::save2mod requires one input argument.');
      end

      % Keep the manager strongly typed so relay restore and other callers
      % cannot silently write malformed structs into cache.
      candidate = varargin{1};
      if ~(isa(candidate, 'pair_symmetry.base.pair_symm_m'))
        error('pair_symmetry::save2mod requires input as a pair_symmetry.base.pair_symm_m object.');
      end

      pair_symm_data = candidate;
      out = int32(0);

    case 'free'
      % Used by testing / stage-restore workflows.
      pair_symm_data = [];
      out = int32(0);

    otherwise
      error('Unknown command: %s', cmd);
  end
end
