function mapping = setsymmGmap(mapping, GWinfor)
  gvec = GWinfor.gvec;
  bz_samp = GWinfor.bz_samp;
  syms = GWinfor.symminfo;

  nrot = syms.nrot;
  ng = gvec.ng;

  indrot = bz_samp.kbz2kibz_ind_rotation;
  indrot = unique(indrot);
  
  g_rot = cell(nrot, 1);
  components = gvec.components;
  % bvec = GWinfor.bvec;



  for id = 1:length(indrot)
    irot = indrot(id);
    mtrx = syms.mtrx{irot};
    indlist = zeros(ng, 1);
    for ig = 1:ng
      gvec_rot = components(ig, :) * mtrx;
      ind = find_vec_in_list(gvec_rot, components);
      if (ind <= 0)
        error('Error in setsymmGmap: gvec_rot not found in components');
      end
      indlist(ig) = ind;
    end
    g_rot{irot} = indlist; 
  end
  mapping.g_rot = g_rot;
  

end % EOF 
