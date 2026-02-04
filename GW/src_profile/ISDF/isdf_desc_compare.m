% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/02/02 ZZ
function flaglist = isdf_desc_compare(isdf_desc1, isdf_desc2)
  %
  flaglist = [true, true, true];
  flagsys = isdf_desc1.desc_sys.equals(isdf_desc2.desc_sys);
  if ~flagsys; return; end
  flagparam = isdf_desc1.desc_param.equals(isdf_desc2.desc_param);
  if ~flagparam; return; end
  % else, compare each 
  ratio1 = isdf_desc1.desc_type1.get('isdf_ratio_type1');
  ratio2 = isdf_desc2.desc_type1.get('isdf_ratio_type1');
  if (ratio1 - ratio2) > 1e-6; flag1 = false; else; flag1 = true; end
  %
  ratio1 = isdf_desc1.desc_type2.get('isdf_ratio_type2');
  ratio2 = isdf_desc2.desc_type2.get('isdf_ratio_type2');
  if (ratio1 - ratio2) > 1e-6; flag2 = false; else; flag2 = true; end
  %
  ratio1 = isdf_desc1.desc_type3.get('isdf_ratio_type3');
  ratio2 = isdf_desc2.desc_type3.get('isdf_ratio_type3');
  if (ratio1 - ratio2) > 1e-6; flag3 = false; else; flag3 = true; end
  %
  flaglist = [~flag1, ~flag2, ~flag3];
  %
  return
end % function