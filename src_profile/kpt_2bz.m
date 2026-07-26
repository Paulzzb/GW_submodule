function kset_inbz = kpt_2bz(kset, fftgrid)
  
  TOL = 1e-6;
  flaghalf = ( abs(round(fftgrid/2) - fftgrid/2) < TOL);
  halffft = round(fftgrid / 2);
  Npoint = size(kset, 1);

  kset_inbz = kset;
  for i = 1:Npoint
    kp = kset(i,:);
    mask = flaghalf & (kp + TOL > halffft);
    kset_inbz(i, mask) = kset_inbz(i, mask) - fftgrid(mask);
  end

end