% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function out = manager(cmd, varargin)
  persistent coulomb_m  

  if nargin==0
    cmd = 'get';
  end

  switch lower(cmd)
    case 'get'
      if isempty(coulomb_m)
        error('coulomb_m not initialized.');
      end
      out = coulomb_m;
      return

    case 'save2mod'
      input = varargin{1};
      if ~isa(input, 'coulomb_m')
        error('coulomb::save2mod required input as a coulomb_m object');
      end
      coulomb_m = input;
      out = 0;

    case 'free'
      coulomb_m = [];

    otherwise
      error('Unknown command')
  end
end