%
% This script shows how to use the scf function to compute the ground state 
% of the bulk Cu using a small unit cell and k-point sampling
%
kssolvpptype('default')
%
%  construct bulk Cu (crystal)
%
%  1. construct atoms
a1 = Atom('Cu');
atomlist = [a1];

%  2. set up the primitive cell
C = [-3.411 0.000 3.411
      0.000 3.411 3.411
     -3.411 3.411 0.000];

%  3. define the coordinates the atoms
xyzlist = [ 0 0 0 ]*C;

% 4. Configure the molecule (crystal)
if (1)
  % Automatic k-points using Monkhorst-Pack grid
  autokpts = [ 4 4 4 0 0 0 ];
  cry = Crystal('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecutwfc',30, 'name','Cu', 'autokpts', autokpts, 'temperature', 300 );
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
    'ecut',30, 'name','Cu', 'temperature', 300, 'kpts', kpts, ...
    'wks', wks);
end

%
% set the options by default and modify some of the options to be passed to the scf function 
opt = setksopt;
opt.maxscfiter = 30;

%
% call the scf function to perform the SCF calcuation
[cry,H,X,info] = scf(cry,opt);

% The following lines are also can be used to compute the band structure of 
% Cu after scf calculation has been performed. They are contain in the 
% plotband function called in the script cu_nscf_k.m
% 
% crykpts = set(cry,'kpts',[0 0 0.25],'wks',[1]); crykpts = finalize(crykpts);
% opt.rho0 = H.rho;
% opt.X0 = X;
% [crykpts,Xkpts,infokpts] = nscf4c(crykpts,opt);
