% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function out = manager(type, cmd, varargin)
  persistent r_lattice d_lattice k q 

  if nargin < 2
    cmd = 'get';
  end

  % Route to the appropriate persistent variable based on type
  switch lower(type)
    case 'r_lat'
      lattice_ptr = r_lattice;
      class_name = 'R_lattice_m';
      
    case 'd_lat'
      lattice_ptr = d_lattice;
      class_name = 'D_lattice_m';
      
    case 'k'
      lattice_ptr = k;
      class_name = 'bz_samp_m';
      
    case 'q'
      lattice_ptr = q;
      class_name = 'bz_samp_m';
      
    otherwise
      error('lattice::manager Unknown type: %s. Use r_lat, d_lat, k, or q', type);
  end

  % Execute command
  switch lower(cmd)
    case 'get'
      if isempty(lattice_ptr)
        error('lattice::%s not initialized. Call lattice.driver first.', type);
      end
      out = lattice_ptr;
      return

    case 'save2mod'     
      % Save to the appropriate persistent variable
      switch lower(type)
        case 'r_lat'
          r_lattice = varargin{1};
        case 'd_lat'
          d_lattice = varargin{1};
        case 'k'
          k = varargin{1};
        case 'q'
          q = varargin{1};
      end
      out = 0;

    case 'free'
      % Clear the appropriate persistent variable
      switch lower(type)
        case 'r_lat'
          r_lattice = [];
        case 'd_lat'
          d_lattice = [];
        case 'k'
          k = [];
        case 'q'
          q = [];
      end
      out = 0;

    otherwise
      error('lattice::manager Unknown command: %s. Use get, save2mod, or free', cmd);
  end
end