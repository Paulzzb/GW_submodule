function v = bandlimV(mol,v)
%
% Usage: v = bandlimV(mol,v);
%
% Purpose:
%   Truncate the potential to be bandlimited.
%
% Input: 
%     mol   a Molecular object
%       v   a 3D array that contains the potential
%
% Output:
%       v   a 3D array that contains the truncated potential
%

idxnz2  = mol.gridrho.idxnz;
vfft = fft3(v);
v = zeros(size(v));
v(idxnz2) = vfft(idxnz2);
v = real(ifft3(v));

end