function mapping = setGomap(mapping, GWinfor)
% function gv = setGomap(gv, nGo)
    % setGomap: generate a mapping, from 
    % nGo          
    % G0_mapping   (:,:) {mustBeInteger}% [nGo x 1]       
  gv = GWinfor.gvec;
  bz_samp = GWinfor.bz_samp;
  iGolist = bz_samp.iGolist;
  iGolist = unique(iGolist);

  nGo = max(iGolist);
  % mapping.nGo = nGo;

  ng = gv.ng;
  n1 = gv.fftgrid(1);
  n2 = gv.fftgrid(2);
  n3 = gv.fftgrid(3);



  Gomapping = cell(nGo,1);
  for id = 1:length(iGolist)
    iGo = iGolist(id);
    tmp = zeros(ng, 1);
    for iG = 1:ng
      G_Go = gv.components(iG,1:3) - gv.components(iGo,1:3);
      G_Go = round(G_Go, 3);
      % G_Go = back2fftzone(G_Go, gvec.fftgrid);
      ind = mod(G_Go(1), n1) + mod(G_Go(2), n2) * n1 ...
            + mod(G_Go(3), n3) * n1 * n2 + 1;
      tmp(iG) = ind;
    end
    Gomapping{iGo} = tmp;
  end

  mapping.Gomapping = Gomapping;
end % EOF
