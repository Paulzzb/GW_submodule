function [aqstemp] = mtxel_sigma(nn, GWinfo, options, sum_range)
% Same as Sigma/mtxel.f90
% Calculate <n| e^{igr}>
% Inputs: 
% nn: 
%
%
%
%
%
%
%
% Output: 
% aqstemp: same as in setGWinfo

% global nv nc  ...
%        ng ...
%        ...
%        ...
%        qk; 

nv = options.nv;
nc = options.nc;

ng = GWinfo.gvec.ng;
qk = GWinfo.qk;
vol = GWinfo.vol;

if nargin < 4
	sum_range = (1 : (nv + nc));
end
% Current version, no gindex needed
gindex            = 1:ng;
Z                 = GWinfo.Z;
im                = sqrt(-1);
gvec              = GWinfo.gvec;
nr = gvec.nfftgridpts; 
% for i = 1:size(Z, 2);
%   Z(:, i) = Z(:, i) / norm(Z(:, i)) * 1;
% end


aqstemp = complex(0.0, 0.0) * zeros(ng, length(sum_range));
fftbox1 = zeros(gvec.fftgrid);
fftbox2 = zeros(gvec.fftgrid);
scale = 1 ./ prod(gvec.fftgrid);

fftbox1 = put_into_fftbox(Z(:, nn), gvec.idxnz, gvec.fftgrid);
% dlmwrite('fftbox1.csv', fftbox1(:));
fftbox1 = nr ./ GWinfo.vol * do_FFT(fftbox1, gvec.fftgrid, 1);
% dlmwrite('fftbox1.csv', fftbox1(:));
fftbox1 = conj(fftbox1);
% dlmwrite('fftbox1.csv', fftbox1(:));

for ind = 1:length(sum_range)
	ind_2 = sum_range(ind);
  fftbox2 = put_into_fftbox(Z(:, ind_2), gvec.idxnz, gvec.fftgrid);
  % dlmwrite('fftbox2.csv', fftbox2(:));
  fftbox2 = nr ./ GWinfo.vol * do_FFT(fftbox2, gvec.fftgrid, 1);
  % dlmwrite('fftbox2.csv', fftbox2(:));
  fftbox2 = fftbox1 .* fftbox2;
	fftbox2 = GWinfo.vol * do_FFT(fftbox2, gvec.fftgrid, 1);
	aqstemp(:, ind) = get_from_fftbox(gvec.idxnz, fftbox2, gvec.fftgrid);
end


return;% main function
end % main function
