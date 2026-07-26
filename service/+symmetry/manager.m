% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function out = manager(cmd, varargin)
  persistent symm_m  

  if nargin==0
    cmd = 'get';
  end

  switch lower(cmd)
    case 'get'
      if isempty(symm_m)
        error('FFT not initialized. Call FFT.init first.');
      end
      out = symm_m;
      return

    case 'save2mod'
      symm_m = varargin{1};
      if isa(symm_m, 'symm_m')
        error('symmetry::save2mod required input as a symm_m object');
      end
      out = 0;

    case 'free'
      symm_m = [];

    otherwise
      error('Unknown command')
  end
end