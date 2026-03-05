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
  flaglist = [true, true, true, true];
  flagsys = isdf_desc1.desc_sys.equals(isdf_desc2.desc_sys);
  if ~flagsys; return; end
  flagparam = isdf_desc1.desc_param.equals(isdf_desc2.desc_param);
  if ~flagparam; return; end
  flag1 = isdf_desc1.desc_type1.equals(isdf_desc2.desc_type1);
  flag2 = isdf_desc1.desc_type2.equals(isdf_desc2.desc_type2);
  flag3 = isdf_desc1.desc_type3.equals(isdf_desc2.desc_type3);
  
  % else, compare each 
  % ratio1 = isdf_desc1.desc_type1.get('isdf_ratio_type1');
  % ratio2 = isdf_desc2.desc_type1.get('isdf_ratio_type1');
  % mlist1 = isdf_desc1.desc_type1.get('mlist');
  % mlist2 = isdf_desc2.desc_type1.get('mlist');
  % nlist1 = isdf_desc1.desc_type1.get('nlist');
  % nlist2 = isdf_desc2.desc_type1.get('nlist');
  % flag1 = true;
  % if (ratio1 - ratio2) > 1e-6;
  %   flag1 = false;
  % elseif all(mlist1 == mlist2)
  %   flag1 = false;
  % elseif all(nlist1 == nlist2)
  %   flag1 = false;
  % end
  %
  % ratio1 = isdf_desc1.desc_type1.get('isdf_ratio_type1');
  % ratio2 = isdf_desc2.desc_type1.get('isdf_ratio_type1');
  % mlist1 = isdf_desc1.desc_type1.get('mlist');
  % mlist2 = isdf_desc2.desc_type1.get('mlist');
  % nlist1 = isdf_desc1.desc_type1.get('nlist');
  % nlist2 = isdf_desc2.desc_type1.get('nlist');
  % flag1 = true;
  % if (ratio1 - ratio2) > 1e-6;
  %   flag1 = false;
  % elseif all(mlist1 == mlist2)
  %   flag1 = false;
  % elseif all(nlist1 == nlist2)
  %   flag1 = false;
  % end
  % ratio1 = isdf_desc1.desc_type2.get('isdf_ratio_type2');
  % ratio2 = isdf_desc2.desc_type2.get('isdf_ratio_type2');
  % if (ratio1 - ratio2) > 1e-6; flag2 = false; else; flag2 = true; end
  % %
  % ratio1 = isdf_desc1.desc_type3.get('isdf_ratio_type3');
  % ratio2 = isdf_desc2.desc_type3.get('isdf_ratio_type3');
  % if (ratio1 - ratio2) > 1e-6; flag3 = false; else; flag3 = true; end
  %
  % -----------------------------------------------------------------
  % Compare vcVnn
  %
  vcmlist1 = isdf_desc1.desc_type1.get('mlist');
  vcmlist2 = isdf_desc2.desc_type1.get('mlist');
  vcnlist1 = isdf_desc1.desc_type1.get('nlist');
  vcnlist2 = isdf_desc2.desc_type1.get('nlist');
  nnmlist1 = isdf_desc1.desc_type3.get('mlist');
  nnmlist2 = isdf_desc2.desc_type3.get('mlist');
  nnnlist1 = isdf_desc1.desc_type3.get('nlist');
  nnnlist2 = isdf_desc2.desc_type3.get('nlist');
  Nvc1 = isdf_desc1.desc_type1.get('Nisdf');
  Nvc2 = isdf_desc2.desc_type1.get('Nisdf');
  Nnn1 = isdf_desc1.desc_type3.get('Nisdf');
  Nnn2 = isdf_desc2.desc_type3.get('Nisdf');
  %
  flag4 = true;
  if ~isequal(vcmlist1, vcmlist2)
    flag4 = false;
  elseif ~isequal(nnmlist1, nnmlist2)
    flag4 = false;
  elseif (Nvc1 ~= Nvc2)
    flag4 = false;
  elseif (Nnn1 ~= Nnn2)
    flag4 = false;
  end

  flaglist = [~flag1, ~flag2, ~flag3, ~flag4];
  %
  return
end % function