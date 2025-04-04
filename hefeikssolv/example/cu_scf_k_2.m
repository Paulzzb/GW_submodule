%
%  construct bulk Cu (crystal)
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

% 4. Configure the molecule (crystal)
if(1)
  % Automatic k-points using Monkhorst-Pack grid
  autokpts = [ 4 4 4 0 0 0 ];
  cry = Crystal('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',50, 'name','Cu', 'autokpts', autokpts, 'temperature', 300 );
else
  % equivalent k-points and weights set manually
  kpts = [
  0.0000  0.0000   0.0000
  0.0000  0.2500   0.0000
  0.0000  0.5000   0.0000
  0.0000  0.7500   0.0000
  0.2500  0.0000   0.0000
  0.2500  0.2500   0.0000
  0.2500  0.5000   0.0000
  0.2500  0.7500   0.0000
  0.5000  0.0000   0.0000
  0.5000  0.2500   0.0000
  0.5000  0.5000   0.0000
  0.5000  0.7500   0.0000
  0.7500  0.0000   0.0000
  0.7500  0.2500   0.0000
  0.7500  0.5000   0.0000
  0.7500  0.7500   0.0000
  0.0000  0.0000   0.2500
  0.0000  0.2500   0.2500
  0.0000  0.5000   0.2500
  0.0000  0.7500   0.2500
  0.2500  0.0000   0.2500
  0.2500  0.2500   0.2500
  0.2500  0.5000   0.2500
  0.2500  0.7500   0.2500
  0.5000  0.0000   0.2500
  0.5000  0.2500   0.2500
  0.5000  0.5000   0.2500
  0.5000  0.7500   0.2500
  0.7500  0.0000   0.2500
  0.7500  0.2500   0.2500
  0.7500  0.5000   0.2500
  0.7500  0.7500   0.2500
  0.0000  0.0000   0.5000
  0.0000  0.2500   0.5000
  0.0000  0.5000   0.5000
  0.0000  0.7500   0.5000
  0.2500  0.0000   0.5000
  0.2500  0.2500   0.5000
  0.2500  0.5000   0.5000
  0.2500  0.7500   0.5000
  0.5000  0.0000   0.5000
  0.5000  0.2500   0.5000
  0.5000  0.5000   0.5000
  0.5000  0.7500   0.5000
  0.7500  0.0000   0.5000
  0.7500  0.2500   0.5000
  0.7500  0.5000   0.5000
  0.7500  0.7500   0.5000
  0.0000  0.0000   0.7500
  0.0000  0.2500   0.7500
  0.0000  0.5000   0.7500
  0.0000  0.7500   0.7500
  0.2500  0.0000   0.7500
  0.2500  0.2500   0.7500
  0.2500  0.5000   0.7500
  0.2500  0.7500   0.7500
  0.5000  0.0000   0.7500
  0.5000  0.2500   0.7500
  0.5000  0.5000   0.7500
  0.5000  0.7500   0.7500
  0.7500  0.0000   0.7500
  0.7500  0.2500   0.7500
  0.7500  0.5000   0.7500
  0.7500  0.7500   0.7500
  ];

  wks = 0.0156250 * ones(size(kpts,1),1);
  cry = Crystal('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',50, 'name','Cu', 'temperature', 300, 'kpts', kpts, ...
    'wks', wks);
end

opt = setksopt;
opt.maxscfiter = 30;

[cry,H,X,info] = scf(cry,opt);

% NSCF calculation
% crykpts = set(cry,'kpts',[0 0 0.25],'wks',[1]); crykpts = finalize(crykpts);
% opt.rho0 = H.rho;
% opt.X0 = X;
% [crykpts,Xkpts,infokpts] = nscf4c(crykpts,opt);
