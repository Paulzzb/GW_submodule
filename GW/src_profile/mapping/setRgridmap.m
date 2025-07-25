function mapping = setRgridmap(mapping, GWinfor)
  bz_samp = GWinfor.bz_samp;
  syms = GWinfor.symminfo;
  fftgridsize = GWinfor.gvec.fftgrid;
  n1 = fftgridsize(1);
  n2 = fftgridsize(2);
  n3 = fftgridsize(3);
  n123 = n1*n2*n3;



  nrot = syms.nrot;
  indrot = bz_samp.kbz2kibz_ind_rotation;
  indrot = unique(indrot);
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main part
  fft_r_rot = cell(nrot, 1);
  for id = 1:length(indrot)
    irot = indrot(id);
    fft_r_rot{irot} = zeros(n123, 1);
  end

  % [gkxind, gkyind, gkzind] = ...
  % ndgrid((0:n1-1)-((0:n1-1) >= n1/2)*n1, ...
  %   (0:n2-1)-((0:n2-1) >= n2/2)*n2, ...
  %   (0:n3-1)-((0:n3-1) >= n3/2)*n3);
  % gkxind = gkxind(:);
  % gkyind = gkyind(:);
  % gkzind = gkzind(:);
  % idlist = 1:length(gkxind);


  for id = 1:length(indrot)
    irot = indrot(id);
    mtrx = syms.mtrx{irot};
    tmpind = zeros(n123, 1);
    count = 0;
    for i3 = 1:n3
      for i2 = 1:n2
        for i1 = 1:n1
          iv = [i1, i2, i3] * mtrx;
          count = count+1;
          i5 = 1 + mod(iv(1), n1) + mod(iv(2), n2)*n1 + mod(iv(3), n3)*n1*n2;
          tmpind(i5) = count;
        end % for i3
      end % for i2
    end % for i1
    fft_r_rot{irot} = tmpind;
  end % for irot
 
  mapping.fft_r_rot = fft_r_rot;

end % EOF 
