% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

classdef symm_m
  properties
    nrot(1, 1)  {mustBeInteger} = 0 
    nsym(1, 1)  {mustBeInteger} = 0 
    is_t_rev(1, 1)  {mustBeInteger} = 1 % 0 or 1
    rot_mtrx_RLU_G  (3, 3, :) = double(0)  % 3x3xnsym
    rot_mtrx_RLU_R  (3, 3, :) = double(0)  % 3x3xnsym
    rot_mtrx_Cart   (3, 3, :) = double(0)  % 3x3xnsym
    inv_rot_index (:, 1) {mustBeInteger} = 0 % nsymx1
  end

  methods
    function obj = symm_m(nsym, nrot, is_t_rev)
    % function obj = symm_m(ob)
      nsym = nsym * (1+is_t_rev);
      obj.nsym = int32( nsym );
      obj.nrot = int32( nrot );
      obj.rot_mtrx_RLU_G = double( zeros(3, 3, nsym) ); 
      obj.rot_mtrx_RLU_R = double( zeros(3, 3, nsym) ); 
      obj.rot_mtrx_Cart  = double( zeros(3, 3, nsym) ); 
      obj.inv_rot_index = int32( zeros(nsym, 1) );
    end % end symm_m 
  end % end method
end % classdef symm_m