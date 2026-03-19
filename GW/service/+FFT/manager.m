% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function out = manager(cmd, varargin)
  persistent data

  if nargin == 0
    cmd = 'get';
  end

  switch lower(cmd)
    case 'get'
      if isempty(data)
        error('FFT not initialized. Call FFT.driver first.');
      end
      out = data;
      return

    case 'save2mod'
      data = varargin{1};
      if ~(isa(data, 'FFT.FFT_m') || isa(data, 'FFT_m'))
        error('FFT::save2mod requires input as a FFT_m object');
      end
      out = 0;

    case 'free'
      data = [];
      out = 0;

    otherwise
      error('Unknown command');
  end
end