function rhog = SCATTER_Bamp(param)
% mtxel_sigma -> Calculate <nn | exp(iGr) | mm>, where mm in sum_range
%     use GWinfo.psir, GWinfo.gvec, and GWinfo.vol
% if sum_range is not specified, sum_range = size(GWinfo.psir, 2)
persistent firsttime DL_vol fftgrid n_spinor

if nargin == 1 && (ischar(param) || (isstring(param) && isscalar(param)))
  cmd = lower(string(param));
  if cmd == "reset"
    firsttime = [];
    DL_vol = [];
    fftgrid = [];
    n_spinor = [];
    rhog = [];
    return
  end
end

if isempty(firsttime)
  firsttime = true;
end

if firsttime
  d_lat_data = lattice.manager('d_lat', 'get');
  DL_vol = d_lat_data.DL_vol;
  fft_data = FFT.manager('get');
  fftgrid = fft_data.fftgrid;
  %
  WF = wave_functions.get();
  n_spinor = WF.n_spinor;
  %
  firsttime = false;
end



iGo = param.qs(1);
irot = param.qs(3);


fft_data = FFT.manager('get');
if irot ~= 1
  r_lat_data = lattice.manager('r_lat', 'get');
end


wf_left = double( wave_functions.WF_apply_symm(param.is) );
wf_right = double( wave_functions.WF_apply_symm(param.os) );



fftbox = conj(wf_left) .* wf_right;



if (n_spinor == 2)
  msg = fprintf("Support for sop calculation is developing\n");
  error(msg);
end

fftbox = reshape(fftbox, fftgrid);
fftbox = DL_vol * do_FFT(fftbox, fftgrid, 1);


if (irot == 1) % no rotation
  rhog = fftbox(  fft_data.G_table(:, iGo)  );
else
  rhog = fftbox(  fft_data.G_table( (r_lat_data.G_rot(:, irot)), iGo)  );
end

rhog = single( rhog );

end % EOF 
