kssolvpptype('default')

%  1. construct atoms
a1 = Atom('Si');
atomlist = [a1 a1 a1 a1 a1 a1 a1 a1];

%  2. set up the primitive cell
C = 10.216*eye(3);

%  3. define the coordinates the atoms
xyzlist = [
    0.000000000         0.000000000         0.000000000
    0.000000000         0.500000000         0.500000000
    0.500000000         0.000000000         0.500000000
    0.500000000         0.500000000         0.000000000
    0.750000000         0.250000000         0.750000000
    0.250000000         0.250000000         0.250000000
    0.250000000         0.750000000         0.750000000
    0.750000000         0.750000000         0.250000000
    ]*C;

%  4. Band structure calculation, from Gamma to X points
kptSta = [0.00  0.00  0.00];  % G
kptEnd = [0.50  0.00  0.00];  % 
nkpts = 21;
kline = (0:nkpts-1)'/(nkpts-1);
kpts = repmat(kptSta,nkpts,1) + kline *(kptEnd-kptSta);

wks = ones(nkpts,1);

crykpts = set(cry,'kpts',kpts,'wks',wks);
crykpts = finalize(crykpts);
opt.rho0 = H.rho;
opt.X0 = X;
[crykpts,Xkpts,infokpts] = nscf4c(crykpts,opt);

% Plot bands (assuming all k-points have the same number of bands)
nbands = length(opt.X0{1}.occ);
eigKpts = reshape(infokpts.Eigvals, nbands, nkpts);

figure
clf
nbnd = 16;
efermi = max(max(eigKpts(1:nbnd,:)));
hold on
for l = 1 : 16
  plot(kline, (eigKpts(l,:)-efermi)*27.211,'-o')
end
hold off


