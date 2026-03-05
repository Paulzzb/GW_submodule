function divide_conquer(GWinfo, config)

  %
  Ngroups = Ne;
  Groups = dad_partition();
  Selected = cell(Ngroups, 1);
  for iGroups = 1:NGroups
    Selected{iGroups} = Select_from_group();
  end
  Subset = concat_Selected();
  Weight4Subset = Weight4Subset();
  %
  

    

end % function