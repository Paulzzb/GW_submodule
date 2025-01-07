% test the element of FFT.
% We try to construct F = KSFFT(mol) explicitly.

cd ../../../
KSSOLV_startup
cd silicon
si32
cd ../src/GW


F_KSSOLV = KSFFT(mol);
vol = mol.vol;
n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;
n123 = mol.n1 * mol.n2 * mol.n3;
nr = n123;
options = GWOptions(mol, []);
% ksinfor = ksinfo(mol, options.Groundstate);
gvec = gvec(mol, 'ecut', mol.ecut);
ng = size(F_KSSOLV);
ii = sqrt(-1);
FMatrix = zeros(nr, ng);

% rindex = ones(3, n1, n2, n3);
% for i = 1:n1
%   rindex(1, i, :, :) = i-1;
% 	rindex(2, :, i, :) = i-1;
% 	rindex(3, :, :, i) = i-1;
% end

% rindex = reshape(rindex, 3, nr);
% rindex = rindex.';

tmparray = zeros(nr, 1);
% for icol = 1:ng
% 	gindex = gvec.components(icol, :) ./ [n1, n2, n3];
%   tmparray = exp(2 * pi * ii * (rindex * gindex'));
%   FMatrix(:, icol) = tmparray; 
% end

FMatrix = FMatrix / vol;

testvecg = zeros(ng, 1);
tol = 1e-10;
tic
for itest = 1:ng
  if itest ~= 1
		testvecg(itest-1) = 0;
	end
	testvecg(itest) = 1;
%	vec1 = FMatrix * testvecg * vol; 
	vec2 = F_KSSOLV' * testvecg * vol;
%	if norm(vec1 - vec2) >= tol
%		fprintf('Not right for itest = %d.\n', itest);
%		error()
%	end
end
toc
tic
for itest = 1:ng
	if itest ~= 1
		testvecg(itest-1) = 0;
	end
	testvecg(itest) = 1;
	fftbox = put_into_fftbox(testvecg, gvec.idxnz, gvec.fftgrid);
	fftbox = do_FFT(fftbox, gvec.fftgrid, 1);
	out = reshape(fftbox, [], 1);
end
toc

