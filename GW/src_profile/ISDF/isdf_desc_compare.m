% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/30 ZZ
function flaglist = isdf_desc_compare(isdf_desc1, isdf_desc2)
  %
  flaglist = [true, true, true];
  flagsys = isdf_desc1.desc_sys.equals(isdf_desc2.desc_sys);
  if ~flagsys; return; end
  flagparam = isdf_desc1.desc_param.equals(isdf_desc2.desc_param);
  if ~flagparam; return; end
  % else, compare each 
  flag1 = isdf_desc1.desc_type1.equals(isdf_desc2.desc_type1);
  flag2 = isdf_desc1.desc_type2.equals(isdf_desc2.desc_type2);
  flag3 = isdf_desc1.desc_type3.equals(isdf_desc2.desc_type3);
  flaglist = [~flag1, ~flag2, ~flag3];
  %
  return
end % function