% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/28 ZZ
function flaglist = isdf_desc_compare(isdf_desc1, isdf_desc2)
  flaglist = [false, false, false];
  flagsys = isdf_desc1.desc_sys.equals(isdf_desc2.desc_sys);
  flagparam = isdf_desc1.desc_param.equals(isdf_desc2.desc_param);
  % if system is not the same, return with all .false.
  if all(flagparam, flagsys)
    return
  end
  % else, compare each 
  flagA = isdf_desc1.descA.equals(isdf_desc2.descA);
  flagB = isdf_desc1.descB.equals(isdf_desc2.descB);
  flagC = isdf_desc1.descC.equals(isdf_desc2.descC);
  flaglist = [flagA, flagB, flagC];
  return
end % function