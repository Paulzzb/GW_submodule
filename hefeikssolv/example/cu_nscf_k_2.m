%
%  construct bulk Cu (crystal). Perform SCF calculations and the band
%  struc
%
kssolvpptype('default')

%  1. construct atoms
a1 = Atom('Cu');
atomlist = [a1];

%  2. set up the primitive cell
%  This matches the definition of FCC cell in
%  http://lampx.tugraz.at/~hadley/ss1/bzones/fcc.php
%  lattice constant from Wannier 90
alat = 3.411 * 2;
C = [ 0.5   0.0   0.5
      0.5   0.5   0.0
      0.0   0.5   0.5 ] * alat;

%  3. define the coordinates the atoms
xyzlist = [
    0.000000000         0.000000000         0.000000000
    ]*C;

%  4. Band structure calculation, from Gamma to X points
kG     = [0.00  0.00  0.00];
kX     = [0.00  0.50  0.50];
kW     = [0.25  0.75  0.50];
kL     = [0.50  0.50  0.50];
kK     = [0.375 0.75  0.375];
nkpts = 41;
kline = (0:nkpts-1)'/(nkpts-1);


kpts      = [];
nkptsTot  = 0;

kptSta = kG;
kptEnd = kX;
nkptsTot = nkptsTot + nkpts;
kpts = [kpts; repmat(kptSta,nkpts,1) + kline *(kptEnd-kptSta);];

kptSta = kX;
kptEnd = kW;
nkptsTot = nkptsTot + nkpts;
kpts = [kpts; repmat(kptSta,nkpts,1) + kline *(kptEnd-kptSta);];

kptSta = kW;
kptEnd = kL;
nkptsTot = nkptsTot + nkpts;
kpts = [kpts; repmat(kptSta,nkpts,1) + kline *(kptEnd-kptSta);];

kptSta = kL;
kptEnd = kG;
nkptsTot = nkptsTot + nkpts;
kpts = [kpts; repmat(kptSta,nkpts,1) + kline *(kptEnd-kptSta);];

kptSta = kG;
kptEnd = kK;
nkptsTot = nkptsTot + nkpts;
kpts = [kpts; repmat(kptSta,nkpts,1) + kline *(kptEnd-kptSta);];


wks = ones(nkptsTot,1);

crykpts = set(cry,'kpts',kpts,'wks',wks);
crykpts = finalize(crykpts);
opt.rho0 = H.rho;
% Assume the same number bands for SCF and non-SCF calculation
opt.X0 = X;
[crykpts,Xkpts,infokpts] = nscf4c(crykpts,opt);

% Plot bands (assuming all k-points have the same number of bands)
nbands = length(opt.X0{1}.occ);
eigKpts = reshape(infokpts.Eigvals, nbands, nkptsTot);

figure
clf
hold on
klineP = kline;
yRange = [13.5 19.5];
for g = 1 : round(nkptsTot / nkpts)
  for l = 1 : nbands
    plot(kline, eigKpts(l,((g-1)*nkpts+1):g*nkpts)*27.211,'b-')
  end
  plot([kline(end) kline(end)], yRange, 'k-');
  kline = kline + 1;
end
hold off
ylim(yRange)
