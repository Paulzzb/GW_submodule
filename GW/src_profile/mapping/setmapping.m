function mapping = setmapping(GWinfor)
  mapping =[];
  mapping = setGomap(mapping, GWinfor);
  mapping = setsymmGmap(mapping, GWinfor);
  mapping = setRgridmap(mapping, GWinfor);

end % EOF